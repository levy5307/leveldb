// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#ifndef STORAGE_LEVELDB_DB_VERSION_EDIT_H_
#define STORAGE_LEVELDB_DB_VERSION_EDIT_H_

#include <set>
#include <utility>
#include <vector>

#include "db/dbformat.h"

namespace leveldb {

class VersionSet;

struct FileMetaData {
  FileMetaData() : refs(0), allowed_seeks(1 << 30), file_size(0) {}

  int refs;
  /** compact之前需要seek的次数 */
  int allowed_seeks;  // Seeks allowed until compaction
  /** number越大的文件越新 */
  uint64_t number;
  /** 文件的大小 */
  uint64_t file_size;    // File size in bytes
  /** 文件中的最小key */
  InternalKey smallest;  // Smallest internal key served by table
  /** 文件中的最大key */
  InternalKey largest;   // Largest internal key served by table
};

/**
 *
 * 表示Version之间的变化，相当于delta增量。Version0 + VersionEdit-->Version1
 *
 * compact过程中会有一系列改变当前Version的操作（FileNumber增加，删除input的sstable，增加
 * 输出的sstable等等），为了缩小Version切换的时间点，将这些操作封装成VersionEdit，compact
 * 完成时，将VersionEdit中的操作一次应用到当前Version即可得到最新状态的Version
 **/
class VersionEdit {
 public:
  VersionEdit() { Clear(); }
  ~VersionEdit() = default;

  void Clear();

  void SetComparatorName(const Slice& name) {
    has_comparator_ = true;
    comparator_ = name.ToString();
  }
  void SetLogNumber(uint64_t num) {
    has_log_number_ = true;
    log_number_ = num;
  }
  void SetPrevLogNumber(uint64_t num) {
    has_prev_log_number_ = true;
    prev_log_number_ = num;
  }
  void SetNextFile(uint64_t num) {
    has_next_file_number_ = true;
    next_file_number_ = num;
  }
  void SetLastSequence(SequenceNumber seq) {
    has_last_sequence_ = true;
    last_sequence_ = seq;
  }
  void SetCompactPointer(int level, const InternalKey& key) {
    compact_pointers_.push_back(std::make_pair(level, key));
  }

  // Add the specified file at the specified number.
  // REQUIRES: This version has not been saved (see VersionSet::SaveTo)
  // REQUIRES: "smallest" and "largest" are smallest and largest keys in file
  void AddFile(int level, uint64_t file, uint64_t file_size,
               const InternalKey& smallest, const InternalKey& largest) {
    FileMetaData f;
    f.number = file;
    f.file_size = file_size;
    f.smallest = smallest;
    f.largest = largest;
    new_files_.push_back(std::make_pair(level, f));
  }

  // Delete the specified "file" from the specified "level".
  void DeleteFile(int level, uint64_t file) {
    deleted_files_.insert(std::make_pair(level, file));
  }

  /** 将VersionEdit成员保存到dst中 */
  void EncodeTo(std::string* dst) const;
  /** 从src中解析VersionEdit成员 */
  Status DecodeFrom(const Slice& src);

  /** 根据VersionEdit成员获取debug string */
  std::string DebugString() const;

 private:
  friend class VersionSet;

  typedef std::set<std::pair<int, uint64_t>> DeletedFileSet;

  std::string comparator_;
  /** log文件的file number */
  uint64_t log_number_;
  uint64_t prev_log_number_;
  /** 下一个可用的 FileNumber */
  uint64_t next_file_number_;
  /** 最后用过的 SequnceNumber */
  SequenceNumber last_sequence_;
  bool has_comparator_;
  bool has_log_number_;
  bool has_prev_log_number_;
  bool has_next_file_number_;
  bool has_last_sequence_;

  /** key-level, value-compact point */
  std::vector<std::pair<int, InternalKey>> compact_pointers_;
  DeletedFileSet deleted_files_;
  std::vector<std::pair<int, FileMetaData>> new_files_;
};

}  // namespace leveldb

#endif  // STORAGE_LEVELDB_DB_VERSION_EDIT_H_
