// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.
//
// The representation of a DBImpl consists of a set of Versions.  The
// newest version is called "current".  Older versions may be kept
// around to provide a consistent view to live iterators.
//
// Each Version keeps track of a set of Table files per level.  The
// entire set of versions is maintained in a VersionSet.
//
// Version,VersionSet are thread-compatible, but require external
// synchronization on all accesses.

#ifndef STORAGE_LEVELDB_DB_VERSION_SET_H_
#define STORAGE_LEVELDB_DB_VERSION_SET_H_

#include <map>
#include <set>
#include <vector>

#include "db/dbformat.h"
#include "db/version_edit.h"
#include "port/port.h"
#include "port/thread_annotations.h"

namespace leveldb {

namespace log {
class Writer;
}

class Compaction;
class Iterator;
class MemTable;
class TableBuilder;
class TableCache;
class Version;
class VersionSet;
class WritableFile;
/**
 * ------------------------------------------------------------------------------------
 * manifest file: {dbname_}/MANIFEST-{manifest_file_number_}
 *   为了重启db后可以恢复退出前的状态，需要将db中的状态保存下来，这些状态信息就保存在manifeest 文件中。
 *   当db出现异常时，为了能够尽可能多的恢复，manifest中不会只保存当前的状态，而是将历史的状态都保存下来。
 *   又考虑到每次状态的完全保存需要的空间和耗费的时间会较多，
 *   当前采用的方式是，只在manifest开始保存完整的状态信息（VersionSet::WriteSnapshot（）），接下来只保存每次compact 产生的操作（VesrionEdit），
 *   重启db时，根据开头的起始状态，依次将后续的VersionEdit replay，即可恢复到退出前的状态（Vesrion）。
 *
 * ------------------------------------------------------------------------------------
 * current file: {dbname_}/CURRENT
 *   保存manifest file name
 *
 * ------------------------------------------------------------------------------------
 * VersionSet:
 *   compact过程中会有一系列改变当前Version的操作（FileNumber增加，删除input的sstable，增加输出的sstable……），
 *   为了缩小Version切换的时间点，将这些操作封装成VersionEdit，compact完成时，将VersionEdit中的操作一次应用到当前Version即可得到最新状态的Version。
 *
 * ------------------------------------------------------------------------------------
 * VersionSet::Builder:
 *   将VersionEdit应用的过程封装成VersionSet::Builder. 主要是更新Version::files_[]
 *
 * ------------------------------------------------------------------------------------
 * VersionEdit:
 *   表示Version之间的变化，相当于delta增量。Version0 + VersionEdit-->Version1
 *
 **/

// Return the smallest index i such that files[i]->largest >= key.
// Return files.size() if there is no such file.
// REQUIRES: "files" contains a sorted list of non-overlapping files.
/** 查找满足files[i]->largest >= key的最小下标i的file, 如果没有找到，则返回files.size() */
int FindFile(const InternalKeyComparator& icmp,
             const std::vector<FileMetaData*>& files, const Slice& key);

// Returns true iff some file in "files" overlaps the user key range
// [*smallest,*largest].
// smallest==nullptr represents a key smaller than all keys in the DB.
// largest==nullptr represents a key largest than all keys in the DB.
// REQUIRES: If disjoint_sorted_files, files[] contains disjoint ranges
//           in sorted order.
/**
 * 如果files中的某个文件中的key集合包含在的[*smallest, *largest]中，则返回true
 * 如果smallest==nullptr, 则代表DB中的最小key
 * 如果largest==nullptr, 则代表DB中的最大key
 **/
bool SomeFileOverlapsRange(const InternalKeyComparator& icmp,
                           bool disjoint_sorted_files,
                           const std::vector<FileMetaData*>& files,
                           const Slice* smallest_user_key,
                           const Slice* largest_user_key);

/**
 * 管理某个版本的所有sstable, db中可能有多个version存在，他们通过链表链接起来
 * 为每个version加入引用计数，读以及解除读操作会将引用计数相应+1/-1，
 * 当version的引用计数为0并且不是最新的version时，该version将会从链表中摘除，
 * 同时，version内的sstable就可以删除了（这些废弃的sstable会在下一次campact完成时删除）
 **/
class Version {
 public:
  // Lookup the value for key.  If found, store it in *val and
  // return OK.  Else return a non-OK status.  Fills *stats.
  // REQUIRES: lock is not held
  struct GetStats {
    FileMetaData* seek_file;
    int seek_file_level;
  };

  // Append to *iters a sequence of iterators that will
  // yield the contents of this Version when merged together.
  // REQUIRES: This version has been saved (see VersionSet::SaveTo)
  /** 获取收集当前版本所有文件的迭代器 */
  void AddIterators(const ReadOptions&, std::vector<Iterator*>* iters);

  /** 在当前版本搜索键值 */
  Status Get(const ReadOptions&, const LookupKey& key, std::string* val,
             GetStats* stats);

  // Adds "stats" into the current state.  Returns true if a new
  // compaction may need to be triggered, false otherwise.
  // REQUIRES: lock is held
  bool UpdateStats(const GetStats& stats);

  // Record a sample of bytes read at the specified internal key.
  // Samples are taken approximately once every config::kReadBytesPeriod
  // bytes.  Returns true if a new compaction may need to be triggered.
  // REQUIRES: lock is held
  bool RecordReadSample(Slice key);

  // Reference count management (so Versions do not disappear out from
  // under live iterators)
  /** 引用 & 解引用 */
  void Ref();
  void Unref();

  /** 获取所有与[begin, end]有overlap的文件集合-->inputs */
  void GetOverlappingInputs(
      int level,
      const InternalKey* begin,  // nullptr means before all keys
      const InternalKey* end,    // nullptr means after all keys
      std::vector<FileMetaData*>* inputs);

  // Returns true iff some file in the specified level overlaps
  // some part of [*smallest_user_key,*largest_user_key].
  // smallest_user_key==nullptr represents a key smaller than all the DB's keys.
  // largest_user_key==nullptr represents a key largest than all the DB's keys.
  /** 判断键值范围与文件是否有交集 */
  bool OverlapInLevel(int level, const Slice* smallest_user_key,
                      const Slice* largest_user_key);

  // Return the level at which we should place a new memtable compaction
  // result that covers the range [smallest_user_key,largest_user_key].
  /** 返回memtable campaction result的存放level */
  int PickLevelForMemTableOutput(const Slice& smallest_user_key,
                                 const Slice& largest_user_key);

  int NumFiles(int level) const { return files_[level].size(); }

  // Return a human readable string that describes this version's contents.
  std::string DebugString() const;

 private:
  friend class Compaction;
  friend class VersionSet;

  class LevelFileNumIterator;

  explicit Version(VersionSet* vset)
      : vset_(vset),
        next_(this),
        prev_(this),
        refs_(0),
        file_to_compact_(nullptr),
        file_to_compact_level_(-1),
        compaction_score_(-1),
        compaction_level_(-1) {}

  Version(const Version&) = delete;
  Version& operator=(const Version&) = delete;

  ~Version();

  /** 获取嵌套双层的two level iterator */
  Iterator* NewConcatenatingIterator(const ReadOptions&, int level) const;

  // Call func(arg, level, f) for every file that overlaps user_key in
  // order from newest to oldest.  If an invocation of func returns
  // false, makes no more calls.
  //
  // REQUIRES: user portion of internal_key == user_key.
  /** 对所有与user_key有overlap的文件执行func函数, 如果func执行返回了false，则停止轮询 */
  void ForEachOverlapping(Slice user_key, Slice internal_key, void* arg,
                          bool (*func)(void*, int, FileMetaData*));

  VersionSet* vset_;  // VersionSet to which this Version belongs
  Version* next_;     // Next version in linked list
  Version* prev_;     // Previous version in linked list
  int refs_;          // Number of live refs to this version

  // List of files per level
  std::vector<FileMetaData*> files_[config::kNumLevels];

  // Next file to compact based on seek stats.
  /** 需要compact的文件(allowed_seeks减少到了0) */
  FileMetaData* file_to_compact_;
  /** file_to_compact_文件的level */
  int file_to_compact_level_;

  // Level that should be compacted next and its compaction score.
  // Score < 1 means compaction is not strictly needed.  These fields
  // are initialized by Finalize().
  double compaction_score_;
  int compaction_level_;
};

/**
 * 整个db的当前状态被VersionSet管理着，其中有当前最新的Version以及其他正在服务的Version链表
 **/
class VersionSet {
 public:
  VersionSet(const std::string& dbname, const Options* options,
             TableCache* table_cache, const InternalKeyComparator*);
  VersionSet(const VersionSet&) = delete;
  VersionSet& operator=(const VersionSet&) = delete;

  ~VersionSet();

  // Apply *edit to the current version to form a new descriptor that
  // is both saved to persistent state and installed as the new
  // current version.  Will release *mu while actually writing to the file.
  // REQUIRES: *mu is held on entry.
  // REQUIRES: no other thread concurrently calls LogAndApply()

  /**
   * 以当前Version为基准构造新的Version，VersionSet::Builder将VersionEdit应用在新Version上，
   * 最后将新Version生效成VersionSet::current_。
   **/
  Status LogAndApply(VersionEdit* edit, port::Mutex* mu)
      EXCLUSIVE_LOCKS_REQUIRED(mu);

  // Recover the last saved descriptor from persistent storage.
  Status Recover(bool* save_manifest);

  // Return the current version.
  Version* current() const { return current_; }

  // Return the current manifest file number
  uint64_t ManifestFileNumber() const { return manifest_file_number_; }

  // Allocate and return a new file number
  uint64_t NewFileNumber() { return next_file_number_++; }

  // Arrange to reuse "file_number" unless a newer file number has
  // already been allocated.
  // REQUIRES: "file_number" was returned by a call to NewFileNumber().
  /** 如果file_number是最新的文件number，则重用该file number */
  void ReuseFileNumber(uint64_t file_number) {
    if (next_file_number_ == file_number + 1) {
      next_file_number_ = file_number;
    }
  }

  // Return the number of Table files at the specified level.
  /** 第level层的文件数量 */
  int NumLevelFiles(int level) const;

  // Return the combined file size of all files at the specified level.
  /** 第level层的文件总大小 */
  int64_t NumLevelBytes(int level) const;

  // Return the last sequence number.
  uint64_t LastSequence() const { return last_sequence_; }

  // Set the last sequence number to s.
  void SetLastSequence(uint64_t s) {
    assert(s >= last_sequence_);
    last_sequence_ = s;
  }

  // Mark the specified file number as used.
  void MarkFileNumberUsed(uint64_t number);

  // Return the current log file number.
  uint64_t LogNumber() const { return log_number_; }

  // Return the log file number for the log file that is currently
  // being compacted, or zero if there is no such log file.
  uint64_t PrevLogNumber() const { return prev_log_number_; }

  // Pick level and inputs for a new compaction.
  // Returns nullptr if there is no compaction to be done.
  // Otherwise returns a pointer to a heap-allocated object that
  // describes the compaction.  Caller should delete the result.
  Compaction* PickCompaction();

  // Return a compaction object for compacting the range [begin,end] in
  // the specified level.  Returns nullptr if there is nothing in that
  // level that overlaps the specified range.  Caller should delete
  // the result.
  Compaction* CompactRange(int level, const InternalKey* begin,
                           const InternalKey* end);

  // Return the maximum overlapping data (in bytes) at next level for any
  // file at a level >= 1.
  int64_t MaxNextLevelOverlappingBytes();

  // Create an iterator that reads over the compaction inputs for "*c".
  // The caller should delete the iterator when no longer needed.
  Iterator* MakeInputIterator(Compaction* c);

  // Returns true iff some level needs a compaction.
  bool NeedsCompaction() const {
    Version* v = current_;
    return (v->compaction_score_ >= 1) || (v->file_to_compact_ != nullptr);
  }

  // Add all files listed in any live version to *live.
  // May also mutate some internal state.
  void AddLiveFiles(std::set<uint64_t>* live);

  // Return the approximate offset in the database of the data for
  // "key" as of version "v".
  uint64_t ApproximateOffsetOf(Version* v, const InternalKey& key);

  // Return a human-readable short (single-line) summary of the number
  // of files per level.  Uses *scratch as backing store.
  struct LevelSummaryStorage {
    char buffer[100];
  };
  const char* LevelSummary(LevelSummaryStorage* scratch) const;

 private:
  class Builder;

  friend class Compaction;
  friend class Version;

  bool ReuseManifest(const std::string& dscname, const std::string& dscbase);

  /** 计算Version *v的compaction level及其compaction score */
  void Finalize(Version* v);

  void GetRange(const std::vector<FileMetaData*>& inputs, InternalKey* smallest,
                InternalKey* largest);

  void GetRange2(const std::vector<FileMetaData*>& inputs1,
                 const std::vector<FileMetaData*>& inputs2,
                 InternalKey* smallest, InternalKey* largest);

  void SetupOtherInputs(Compaction* c);

  // Save current contents to *log
  Status WriteSnapshot(log::Writer* log);

  void AppendVersion(Version* v);

  Env* const env_;
  const std::string dbname_;
  const Options* const options_;
  TableCache* const table_cache_;
  const InternalKeyComparator icmp_;
  /** 下一个可用的FileNumber */
  uint64_t next_file_number_;
  /** manifest文件的FileNumber */
  uint64_t manifest_file_number_;
  /** 最后用过的 SequnceNumber */
  uint64_t last_sequence_;
  /** log 文件的 FileNumber */
  uint64_t log_number_;
  uint64_t prev_log_number_;  // 0 or backing store for memtable being compacted

  // Opened lazily
  WritableFile* descriptor_file_;
  log::Writer* descriptor_log_;
  /** version链表的表头（双向循环链表） */
  Version dummy_versions_;  // Head of circular doubly-linked list of versions.
  /** 当前version: 指向头结点的next，也就是最新插入的结点 */
  Version* current_;        // == dummy_versions_.prev_

  // Per-level key at which the next compaction at that level should start.
  // Either an empty string, or a valid InternalKey.
  std::string compact_pointer_[config::kNumLevels];
};

// A Compaction encapsulates information about a compaction.
class Compaction {
 public:
  ~Compaction();

  // Return the level that is being compacted.  Inputs from "level"
  // and "level+1" will be merged to produce a set of "level+1" files.
  int level() const { return level_; }

  // Return the object that holds the edits to the descriptor done
  // by this compaction.
  VersionEdit* edit() { return &edit_; }

  // "which" must be either 0 or 1
  int num_input_files(int which) const { return inputs_[which].size(); }

  // Return the ith input file at "level()+which" ("which" must be 0 or 1).
  FileMetaData* input(int which, int i) const { return inputs_[which][i]; }

  // Maximum size of files to build during this compaction.
  uint64_t MaxOutputFileSize() const { return max_output_file_size_; }

  // Is this a trivial compaction that can be implemented by just
  // moving a single input file to the next level (no merging or splitting)
  /** trivial move: 只需要把一个文件移动到下一个level */
  bool IsTrivialMove() const;

  // Add all inputs to this compaction as delete operations to *edit.
  /** 将所有本次compaction操作的input文件，在edit中标记为删除 */
  void AddInputDeletions(VersionEdit* edit);

  // Returns true if the information we have available guarantees that
  // the compaction is producing data in "level+1" for which no data exists
  // in levels greater than "level+1".
  /** 如果user_key在level_+1以上的的level中不存在，返回true(user_key是compaction产生的在level+1中的) */
  bool IsBaseLevelForKey(const Slice& user_key);

  // Returns true if we should stop building the current output
  // before processing "internal_key".
  /** 在处理internal_key之前，如果当前与grandparent层产生overlap的size超过阈值, 返回true */
  bool ShouldStopBefore(const Slice& internal_key);

  // Release the input version for the compaction, once the compaction
  // is successful.
  /** 当compaction成功执行完后，释放掉input_verison_ */
  void ReleaseInputs();

 private:
  friend class Version;
  friend class VersionSet;

  Compaction(const Options* options, int level);

  /** 要compact的level */
  int level_;
  /** 生成的sstable的最大size */
  uint64_t max_output_file_size_;
  /** compact时当前的version */
  Version* input_version_;
  /** 记录compact过程中的操作 */
  VersionEdit edit_;

  // Each compaction reads inputs from "level_" and "level_+1"
  /**
   * inputs_[0]是level_的sstable文件信息
   * inputs_[1]是level_+1的sstable文件信息
   **/
  std::vector<FileMetaData*> inputs_[2];  // The two sets of inputs

  // State used to check for number of overlapping grandparent files
  // (parent == level_ + 1, grandparent == level_ + 2)
  /**
   * 位于level_+2，并且与compact的key-range有overlap的sstable。
   * 保存grandparents_是因为compact最终会生成一系列level_+1的sstable，
   * 而如果生成的sstable与level_+2中有过多的overlap的话，当compact
   * level_+1时，会产生过多的merge，为了尽量避免这种情况，compact过程中
   * 需要检查与level-n+2中产生overlap的size并与阈值kMaxGrandParentOverlapBytes做比较，
   * 以便提前中止compact。
   **/
  std::vector<FileMetaData*> grandparents_;
  /** 记录compact时grandparents_中已经overlap的index */
  size_t grandparent_index_;  // Index in grandparent_starts_
  bool seen_key_;             // Some output key has been seen
  /**  记录已经overlap的累计size */
  int64_t overlapped_bytes_;  // Bytes of overlap between current output
                              // and grandparent files

  // State for implementing IsBaseLevelForKey

  // level_ptrs_ holds indices into input_version_->levels_: our state
  // is that we are positioned at one of the file ranges for each
  // higher level than the ones involved in this compaction (i.e. for
  // all L >= level_ + 2).
  /**
   * compact时key的遍历是顺序的，所以每次检查从上一次检查结束的地方开始即可，
   * level_ptrs_[i]中就记录了input_version_->levels_[i]中，上一次比较结束的sstable的容器下标。
   * (compact时，当key的ValueType是kTypeDeletion时，要检查其在level-n+1以上是否存在（IsBaseLevelForKey()）来决定是否丢弃掉该key)
   *
   * 这里真正需要的是level>level_+2的层次中的、某一文件的内容下标。
   **/
  size_t level_ptrs_[config::kNumLevels];
};

}  // namespace leveldb

#endif  // STORAGE_LEVELDB_DB_VERSION_SET_H_
