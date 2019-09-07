// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#include "db/db_iter.h"

#include "db/db_impl.h"
#include "db/dbformat.h"
#include "db/filename.h"
#include "leveldb/env.h"
#include "leveldb/iterator.h"
#include "port/port.h"
#include "util/logging.h"
#include "util/mutexlock.h"
#include "util/random.h"

namespace leveldb {

#if 0
static void DumpInternalIter(Iterator* iter) {
  for (iter->SeekToFirst(); iter->Valid(); iter->Next()) {
    ParsedInternalKey k;
    if (!ParseInternalKey(iter->key(), &k)) {
      fprintf(stderr, "Corrupt '%s'\n", EscapeString(iter->key()).c_str());
    } else {
      fprintf(stderr, "@ '%s'\n", k.DebugString().c_str());
    }
  }
}
#endif

namespace {

// Memtables and sstables that make the DB representation contain
// (userkey,seq,type) => uservalue entries.  DBIter
// combines multiple entries for the same userkey found in the DB
// representation into a single entry while accounting for sequence
// numbers, deletion markers, overwrites, etc.
/**
 *  ———————————————————————————————————————————————————
 * | user key | sequence number(7 bytes) | type(1 byte)|
 *  ----------------------------------------------------
 * key由三部分组成：
 *  1.用户定义的key
 *  2.序列号：leveldb中，每次写操作都有一个sequence number，标志写入的先后顺序。由于在leveldb中，可能会有多条相同的
 *           key的数据项同时存储在数据库中，因此需要一个序列号来标识这些序号的新旧情况。序列号最大的数据项为最新值
 *  3.标志本条数据项的类型，是更新还是删除：
 *
 *  db中key的排列顺序，key由小到大，相同的key的情况下，sequence number由大到小
 *  例如: 下图中，key1 < key2, seq1 > seq2
 *  (key1, seq1, type) (key1, seq2, type) (key2, seq1, type) (key2, seq2, type) (key3, seq3, type)
 **/
class DBIter : public Iterator {
 public:
  // Which direction is the iterator currently moving?
  // (1) When moving forward, the internal iterator is positioned at
  //     the exact entry that yields this->key(), this->value()
  // (2) When moving backwards, the internal iterator is positioned
  //     just before all entries whose user key == this->key().
  enum Direction { kForward, kReverse };

  DBIter(DBImpl* db, const Comparator* cmp, Iterator* iter, SequenceNumber s,
         uint32_t seed)
      : db_(db),
        user_comparator_(cmp),
        iter_(iter),
        sequence_(s),
        direction_(kForward),
        valid_(false),
        rnd_(seed),
        bytes_until_read_sampling_(RandomCompactionPeriod()) {}

  DBIter(const DBIter&) = delete;
  DBIter& operator=(const DBIter&) = delete;

  ~DBIter() override { delete iter_; }
  bool Valid() const override { return valid_; }
  /**
   * 返回当前key，如果方向是kFoward，iter_指向的就是当前entry，直接取就可以了；
   * 否则当前key保存在saved_key_中
   **/
  Slice key() const override {
    assert(valid_);
    return (direction_ == kForward) ? ExtractUserKey(iter_->key()) : saved_key_;
  }
    /**
     * 返回当前value，如果方向是kFoward，iter_指向的就是当前entry，直接取就可以了；
     * 否则当前value保存在saved_value_中
     **/
  Slice value() const override {
    assert(valid_);
    return (direction_ == kForward) ? iter_->value() : saved_value_;
  }
  Status status() const override {
    if (status_.ok()) {
      return iter_->status();
    } else {
      return status_;
    }
  }

  void Next() override;
  void Prev() override;
  void Seek(const Slice& target) override;
  void SeekToFirst() override;
  void SeekToLast() override;

 private:
  void FindNextUserEntry(bool skipping, std::string* skip);
  void FindPrevUserEntry();
  bool ParseKey(ParsedInternalKey* key);

  inline void SaveKey(const Slice& k, std::string* dst) {
    dst->assign(k.data(), k.size());
  }

  inline void ClearSavedValue() {
    if (saved_value_.capacity() > 1048576) {
      std::string empty;
      swap(empty, saved_value_);
    } else {
      saved_value_.clear();
    }
  }

  // Picks the number of bytes that can be read until a compaction is scheduled.
  size_t RandomCompactionPeriod() {
    return rnd_.Uniform(2 * config::kReadBytesPeriod);
  }

  DBImpl* db_;
  const Comparator* const user_comparator_;
  Iterator* const iter_;
  SequenceNumber const sequence_;
  Status status_;
  std::string saved_key_;    // == current key when direction_==kReverse
  std::string saved_value_;  // == current raw value when direction_==kReverse
  Direction direction_;
  bool valid_;
  Random rnd_;
  size_t bytes_until_read_sampling_;
};

inline bool DBIter::ParseKey(ParsedInternalKey* ikey) {
  Slice k = iter_->key();

  /**
   * 获取当前iter_所指向的key，如果bytes_untili_read_sampling < bytes_read, 则记录sample,
   * 并随机增加bytes_until_read_sampling的值
   * */
  size_t bytes_read = k.size() + iter_->value().size();
  while (bytes_until_read_sampling_ < bytes_read) {
    bytes_until_read_sampling_ += RandomCompactionPeriod();
    db_->RecordReadSample(k);
  }

  /** 对byte_until_read_sampling减去bytes_read */
  assert(bytes_until_read_sampling_ >= bytes_read);
  bytes_until_read_sampling_ -= bytes_read;

  /** 解析获取ParsedInternalKey类型的ikey */
  if (!ParseInternalKey(k, ikey)) {
    status_ = Status::Corruption("corrupted internal key in DBIter");
    return false;
  } else {
    return true;
  }
}

void DBIter::Next() {
  assert(valid_);

  /** 对于next操作，如果当前direction不是kForward, 则需要先切换方向 */
  if (direction_ == kReverse) {  // Switch directions?
    direction_ = kForward;
    // iter_ is pointing just before the entries for this->key(),
    // so advance into the range of entries for this->key() and then
    // use the normal skipping code below.
    /**
     *  1. 如果与前一次操作 dirction_一致，直接 FindNextUserEntry()
     *  2. 否则，前一次的 Prev（）使 iter_定位在当前key的前一个，先 iter_->Next()（如果已经 Prev（）遍历完了，则iter_->SeekToFirst()）,
     *     回到当前位置，然后再FindNextUserEntry()。
     *  备注：这里的saved_key_保存的是当前key, 当调用FindNextUserEntry时可以跳过所有与当前key相同的entry
     **/
    if (!iter_->Valid()) {
      iter_->SeekToFirst();
    } else {
      iter_->Next();
    }
    if (!iter_->Valid()) {
      valid_ = false;
      saved_key_.clear();
      return;
    }
    // saved_key_ already contains the key to skip past.
  } else {
    // Store in saved_key_ the current key so we skip it below.
    /** 解析当前iter_指向的key，获取其userkey-->saved_key_ */
    SaveKey(ExtractUserKey(iter_->key()), &saved_key_);

    // iter_ is pointing to current key. We can now safely move to the next to
    // avoid checking current key.
    /**
     * 1.并给FindNextUserEntry的skipping设置为true, 并将当前iter_指向的key保存到saved_key中，
     * 2.iter_安全后移
     * 3.将二者传递给FindNextUserEntry，用于跳过与iter_指向的相同的key
     * 备注：这里也可以不执行else里的这些语句，直接调FindNextUserEntry就可以了。
     **/
    iter_->Next();
    if (!iter_->Valid()) {
      valid_ = false;
      saved_key_.clear();
      return;
    }
  }

  FindNextUserEntry(true, &saved_key_);
}

/**
 * 正向遍历，首次遇到的key就是key的最终状态（SequnceNumber更大），处理简单:
 *   1.读取iter_指向的key，如果该key类型是kTypeDeletion, 则说明该key被删除了，则继续跳过所有相同的key，
 *   2.直到获取到下一个不同的key，重新执行1
 * 返回时，iter_指向下一个有效的key（对于相同的key，指向sequence最大的那个）
 **/
void DBIter::FindNextUserEntry(bool skipping, std::string* skip) {
  // Loop until we hit an acceptable entry to yield
  assert(iter_->Valid());
  assert(direction_ == kForward);
  do {
    ParsedInternalKey ikey;
    /** 解析获取当前key --> ikey */
    if (ParseKey(&ikey) && ikey.sequence <= sequence_) {
      switch (ikey.type) {
        /** 说明该key被删除了，保存该key到skip，在后续读时，对于相同的key直接跳过 */
        case kTypeDeletion:
          // Arrange to skip all upcoming entries for this key since
          // they are hidden by this deletion.
          SaveKey(ikey.user_key, skip);
          skipping = true;
          break;
        case kTypeValue:
          /** 该key被删除了，直接跳过 */
          if (skipping &&
              user_comparator_->Compare(ikey.user_key, *skip) <= 0) {
            // Entry hidden
          } else {
            valid_ = true;
            saved_key_.clear();
            return;
          }
          break;
      }
    }
    iter_->Next();
  } while (iter_->Valid());
  saved_key_.clear();
  valid_ = false;
}

void DBIter::Prev() {
  assert(valid_);

  /**
   * 方向不是kBackward，切换方向
   * 由于方向是kForward，说明当前iter_指向当前key，需要回退到当前key的前一个
   **/
  if (direction_ == kForward) {  // Switch directions?
    // iter_ is pointing at the current entry.  Scan backwards until
    // the key changes so we can use the normal reverse scanning code.
    assert(iter_->Valid());  // Otherwise valid_ would have been false
    SaveKey(ExtractUserKey(iter_->key()), &saved_key_);
    while (true) {
      /** 回退一步 */
      iter_->Prev();

      /** 回退到了开头了，则无法找Prev了，返回 */
      if (!iter_->Valid()) {
        valid_ = false;
        saved_key_.clear();
        ClearSavedValue();
        return;
      }

      /** 找到了当前key的前一个，直接返回 */
      if (user_comparator_->Compare(ExtractUserKey(iter_->key()), saved_key_) <
          0) {
        break;
      }
    }
    direction_ = kReverse;
  }

  FindPrevUserEntry();
}

/**
 * 反向遍历时，对于一个key，最后遍历到的才是其最终状态，所以必须遍历到该key的前一个，才能确定该key已经全部处理过，并获得其最终状态。
 * 这时iter_并不位于当前key的位置，所以需要saved_key_/save_value_来保存当前的key/value。
 **/
void DBIter::FindPrevUserEntry() {
  assert(direction_ == kReverse);

  ValueType value_type = kTypeDeletion;
  if (iter_->Valid()) {
    do {
      ParsedInternalKey ikey;
      if (ParseKey(&ikey) && ikey.sequence <= sequence_) {
        /** 找到了比saved_key_（当前key）更小的key，则说明遍历到了当前key的前一个，可以返回了 */
        if ((value_type != kTypeDeletion) &&
            user_comparator_->Compare(ikey.user_key, saved_key_) < 0) {
          // We encountered a non-deleted value in entries for previous keys,
          break;
        }

        /********************************************************
         *** 执行到这里，说明当前遍历的iter的key与saved_key_相同  ***
         ********************************************************/

        /** 如果当前遍历的entry的value_type是kTypeDeletion, 则清除slaved_key_和slaved_value_ */
        value_type = ikey.type;
        if (value_type == kTypeDeletion) {
          saved_key_.clear();
          ClearSavedValue();
        } else {
          /** 由于是向前遍历，前面的entry的sequence num更大，代表数据更新，所以更新saved_key和saved_value */
          Slice raw_value = iter_->value();
          if (saved_value_.capacity() > raw_value.size() + 1048576) {
            std::string empty;
            swap(empty, saved_value_);
          }
          SaveKey(ExtractUserKey(iter_->key()), &saved_key_);
          saved_value_.assign(raw_value.data(), raw_value.size());
        }
      }
      iter_->Prev();
    } while (iter_->Valid());
  }

  if (value_type == kTypeDeletion) {
    /**
     * value_type == kTypeDeletion这里说明是由于while(iter_->Valid())不满足而退出的(因为如果是break退出的，value_type肯定不等于kTypeDeletion)，
     * 说明执行到了最开头，并且当前key已经被删掉了
     **/
    // End
    valid_ = false;
    saved_key_.clear();
    ClearSavedValue();
    direction_ = kForward;
  } else {
    valid_ = true;
  }
}

void DBIter::Seek(const Slice& target) {
  direction_ = kForward;
  ClearSavedValue();
  saved_key_.clear();
  /** 将target保存到saved_key_中 */
  AppendInternalKey(&saved_key_,
                    ParsedInternalKey(target, sequence_, kValueTypeForSeek));

  /** 定位到saved_key_ */
  iter_->Seek(saved_key_);
  if (iter_->Valid()) {
    /** 正向遍历，skipping设置为false表示不跳过saved_key_, 用于跳过被删除的entry(value_type != kTypeDeletion) */
    FindNextUserEntry(false, &saved_key_ /* temporary storage */);
  } else {
    valid_ = false;
  }
}

void DBIter::SeekToFirst() {
  direction_ = kForward;
  ClearSavedValue();
  iter_->SeekToFirst();
  if (iter_->Valid()) {
    /**
     * 正向遍历，skipping设置为false表示不跳过saved_key_,
     * 用于跳过被删除的entry(value_type != kTypeDeletion)、定位到当前key的iter_
     **/
    FindNextUserEntry(false, &saved_key_ /* temporary storage */);
  } else {
    valid_ = false;
  }
}

void DBIter::SeekToLast() {
  direction_ = kReverse;
  ClearSavedValue();
  iter_->SeekToLast();
  /** 越过被删除的entry(value_type != kTypeDeletion)、将iter_定位到last entry的前一个 */
  FindPrevUserEntry();
}

}  // anonymous namespace

Iterator* NewDBIterator(DBImpl* db, const Comparator* user_key_comparator,
                        Iterator* internal_iter, SequenceNumber sequence,
                        uint32_t seed) {
  return new DBIter(db, user_key_comparator, internal_iter, sequence, seed);
}

}  // namespace leveldb
