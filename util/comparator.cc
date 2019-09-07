// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#include "leveldb/comparator.h"

#include <algorithm>
#include <cstdint>
#include <string>
#include <type_traits>

#include "leveldb/slice.h"
#include "util/logging.h"
#include "util/no_destructor.h"

namespace leveldb {

Comparator::~Comparator() = default;

namespace {
class BytewiseComparatorImpl : public Comparator {
 public:
  BytewiseComparatorImpl() = default;

  const char* Name() const override { return "leveldb.BytewiseComparator"; }

  /** 重构Compare函数 */
  int Compare(const Slice& a, const Slice& b) const override {
    return a.compare(b);
  }

  /**
   * 找出[start, limit]区间内的一个最小的短串
   * 例如：start = [80, 86, 88, 99, 20]
   *      limit = [80, 86, 90]
   * 找出的最小短串是[80, 86, 89]
   */
  void FindShortestSeparator(std::string* start,
                             const Slice& limit) const override {
    // Find length of common prefix
    size_t min_length = std::min(start->size(), limit.size());
    size_t diff_index = 0;
    while ((diff_index < min_length) &&
           ((*start)[diff_index] == limit[diff_index])) {
      diff_index++;
    }

    if (diff_index >= min_length) {
      // Do not shorten if one string is a prefix of the other
    } else {
      uint8_t diff_byte = static_cast<uint8_t>((*start)[diff_index]);
      if (diff_byte < static_cast<uint8_t>(0xff) &&
          diff_byte + 1 < static_cast<uint8_t>(limit[diff_index])) {
        (*start)[diff_index]++;
        start->resize(diff_index + 1);
        assert(Compare(*start, limit) < 0);
      }
    }
  }

  /**
   * 寻找第一个可以+1的(也就是 < 255)byte, 然后对该byte加一，并将key剪切到该byte
   * 也就是找到该key的最小的后继，例如 [255, 255，121，80]的最小的后继是[255, 255, 122]
   * */
  void FindShortSuccessor(std::string* key) const override {
    // Find first character that can be incremented
    size_t n = key->size();
    for (size_t i = 0; i < n; i++) {
      const uint8_t byte = (*key)[i];
      if (byte != static_cast<uint8_t>(0xff)) {
        (*key)[i] = byte + 1;
        key->resize(i + 1);
        return;
      }
    }
    // *key is a run of 0xffs.  Leave it alone.
  }
};
}  // namespace

/** 这种singleton的实现方式是线程安全的 */
const Comparator* BytewiseComparator() {
  static NoDestructor<BytewiseComparatorImpl> singleton;
  return singleton.get();
}

}  // namespace leveldb
