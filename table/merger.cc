// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#include "table/merger.h"

#include "leveldb/comparator.h"
#include "leveldb/iterator.h"
#include "table/iterator_wrapper.h"

namespace leveldb {

/**
 * 结构图：
 *  管理n个iterator（child）, 每个iterator指向一列元素。
 *  支持:
 *      1.查找这n个iterator所有指向的元素的最小key
 *      1.查找这n个iterator所有指向的元素的最大key
 *      2.查找当前指向元素的下一个(比其key值大)元素
 *      3.查找当前指向元素的前一个(比其key值小)元素
 *
 *   __  升序         ___                  ___
 *  |__|  |          |__|            .    |__|
 *  |__|  |          |__|            .    |__|<--child(n-1)
 *  |__|  V          |__|            .    |__|
 *  |__| <--child0   |__|            .    |__|
 *  |__|             |__|            .    |__|
 *  |__|             |__|<--child1   .    |__|
 *  |__|             |__|                 |__|
 **/
namespace {
class MergingIterator : public Iterator {
 public:
  MergingIterator(const Comparator* comparator, Iterator** children, int n)
      : comparator_(comparator),
        children_(new IteratorWrapper[n]),
        n_(n),
        current_(nullptr),
        direction_(kForward) {
    for (int i = 0; i < n; i++) {
      children_[i].Set(children[i]);
    }
  }

  ~MergingIterator() override { delete[] children_; }

  bool Valid() const override { return (current_ != nullptr); }

  void SeekToFirst() override {
    /** 所有的child都seek to fist, 并令current_指向最小的child */
    for (int i = 0; i < n_; i++) {
      children_[i].SeekToFirst();
    }
    FindSmallest();

    /** 令行进方向=向前 */
    direction_ = kForward;
  }

  void SeekToLast() override {
    /** 所有的child都seek to last, 并令current_指向最大的child */
    for (int i = 0; i < n_; i++) {
      children_[i].SeekToLast();
    }
    FindLargest();

    /** 令行进方向=向后 */
    direction_ = kReverse;
  }

  void Seek(const Slice& target) override {
    /** 所有的child都指向target, 并令current_指向最小的child */
    for (int i = 0; i < n_; i++) {
      children_[i].Seek(target);
    }
    FindSmallest();

    /** 令行进方向=向前 */
    direction_ = kForward;
  }

  /** 用于在所有child中寻找，比current_的key大的下一个iterator */
  void Next() override {
    assert(Valid());

    // Ensure that all children are positioned after key().
    // If we are moving in the forward direction, it is already
    // true for all of the non-current_ children since current_ is
    // the smallest child and key() == current_->key().  Otherwise,
    // we explicitly position the non-current_ children.
    /** 如果方向是kBackward，则说明上次操作方向是kBackward, 这时候current_指向的元素比所有childs都大 */
    if (direction_ != kForward) {
      for (int i = 0; i < n_; i++) {
        IteratorWrapper* child = &children_[i];
        if (child != current_) {
          child->Seek(key());
          if (child->Valid() &&
              comparator_->Compare(key(), child->key()) == 0) {
            child->Next();
          }
        }
      }
      direction_ = kForward;
    }

    current_->Next();
    FindSmallest();
  }

  void Prev() override {
    assert(Valid());

    // Ensure that all children are positioned before key().
    // If we are moving in the reverse direction, it is already
    // true for all of the non-current_ children since current_ is
    // the largest child and key() == current_->key().  Otherwise,
    // we explicitly position the non-current_ children.
    /**
     *  因为有多个Iterator存在，需要记录前一次做的是何种方向的操作, 判断这一次操作的方向是否和前一次一致，来做不同的处理。
     *  比如，如果之前一步做了Next(), 那么direction_就是kForward, current_定位到的一定是children_中满足条件最小的(顺序由小到大)，
     *  其他的Iterator已经定位到大于当前key的位置(除非Iterator已经end), 这时继续做Next(), 只需要current_->Next()，然后在children_中选出大于当前key且最小的即可.
     *  但如果做Prev(), 其他的 Iterator可能位于大于当前key的位置，所以必须先让所有的Iterator都定位到小于当前key的位置(Ierator中不存在key，就SeekToLast()),
     *  然后选出小于当前key且最大的
     **/
    if (direction_ != kReverse) {
      for (int i = 0; i < n_; i++) {
        IteratorWrapper* child = &children_[i];
        if (child != current_) {
          child->Seek(key());
          if (child->Valid()) {
            // Child is at first entry >= key().  Step back one to be < key()
            child->Prev();
          } else {
            // Child has no entries >= key().  Position at last entry.
            child->SeekToLast();
          }
        }
      }
      direction_ = kReverse;
    }

    current_->Prev();
    FindLargest();
  }

  Slice key() const override {
    assert(Valid());
    return current_->key();
  }

  Slice value() const override {
    assert(Valid());
    return current_->value();
  }

  /** 只要有一个child的status不是ok, 就返回not ok */
  Status status() const override {
    Status status;
    for (int i = 0; i < n_; i++) {
      status = children_[i].status();
      if (!status.ok()) {
        break;
      }
    }
    return status;
  }

 private:
  // Which direction is the iterator moving?
  /**
   *  因为有多个Iterator存在，需要记录前一次做的是何种方向的操作, 判断这一次操作的方向是否和前一次一致，来做不同的处理。
   *  比如，如果做了Next()，current_定位到的一定是children_中满足条件最小的，其他的Iterator已经定位到大于当前key的位置(除非Iterator已经end），
   *  这时继续做Next(), 只需要current_->Next()，然后在children_中选出大于当前key且最小的即可.
   *  但如果做Prev(), 其他的 Iterator可能位于大于当前key的位置，所以必须先让所有的Iterator都定位到小于当前key的位置(Ierator中不存在key，就SeekToLast()),
   *  然后选出小于当前key且最大的
   **/
  enum Direction { kForward, kReverse };

  void FindSmallest();
  void FindLargest();

  const Comparator* comparator_;
  // We might want to use a heap in case there are lots of children.
  // For now we use a simple array since we expect a very small number
  // of children in leveldb.
  IteratorWrapper* children_;
  int n_;
  IteratorWrapper* current_;
  Direction direction_;
};

/** 查找key最小的child, 令curren_指向它 */
void MergingIterator::FindSmallest() {
  IteratorWrapper* smallest = nullptr;
  for (int i = 0; i < n_; i++) {
    IteratorWrapper* child = &children_[i];
    if (child->Valid()) {
      if (smallest == nullptr) {
        smallest = child;
      } else if (comparator_->Compare(child->key(), smallest->key()) < 0) {
        smallest = child;
      }
    }
  }
  current_ = smallest;
}

/** 查找key最大的child, 令curren_指向它 */
void MergingIterator::FindLargest() {
  IteratorWrapper* largest = nullptr;
  for (int i = n_ - 1; i >= 0; i--) {
    IteratorWrapper* child = &children_[i];
    if (child->Valid()) {
      if (largest == nullptr) {
        largest = child;
      } else if (comparator_->Compare(child->key(), largest->key()) > 0) {
        largest = child;
      }
    }
  }
  current_ = largest;
}
}  // namespace

/** 创建一个新的MergingIterator */
Iterator* NewMergingIterator(const Comparator* comparator, Iterator** children,
                             int n) {
  assert(n >= 0);
  if (n == 0) {
    return NewEmptyIterator();
  } else if (n == 1) {
    return children[0];
  } else {
    return new MergingIterator(comparator, children, n);
  }
}

}  // namespace leveldb
