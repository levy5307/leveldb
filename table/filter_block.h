// Copyright (c) 2012 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.
//
// A filter block is stored near the end of a Table file.  It contains
// filters (e.g., bloom filters) for all data blocks in the table combined
// into a single filter block.

#ifndef STORAGE_LEVELDB_TABLE_FILTER_BLOCK_H_
#define STORAGE_LEVELDB_TABLE_FILTER_BLOCK_H_

#include <stddef.h>
#include <stdint.h>

#include <string>
#include <vector>

#include "leveldb/slice.h"
#include "util/hash.h"

namespace leveldb {

class FilterPolicy;

// A FilterBlockBuilder is used to construct all of the filters for a
// particular Table.  It generates a single string which is stored as
// a special block in the Table.
//
// The sequence of calls to FilterBlockBuilder must match the regexp:
//      (StartBlock AddKey*)* Finish
/**
 * filter block就是meta block, 存储了block的filter数据, 用于加快查询的速度
 * filter block格式:
 *    |      filter data 1
 *    |      filter data 2
 *    |           ...
 *    |      filter data n
 *    |     filter offset 1
 *    |     filter offset 2
 *    |           ...
 *    |     filter offset n
 *    |  beginning of filter offset
 *    V          base
 * */
class FilterBlockBuilder {
 public:
  explicit FilterBlockBuilder(const FilterPolicy*);

  FilterBlockBuilder(const FilterBlockBuilder&) = delete;
  FilterBlockBuilder& operator=(const FilterBlockBuilder&) = delete;

  void StartBlock(uint64_t block_offset);
  void AddKey(const Slice& key);
  Slice Finish();

 private:
  void GenerateFilter();

  const FilterPolicy* policy_;
  /** 缓存了所有key */
  std::string keys_;             // Flattened key contents
  /** 缓存所有的key的偏移 */
  std::vector<size_t> start_;    // Starting index in keys_ of each key
  /** 缓存了当前的计算结果，最后通过finish，将result_以及filter_offsets组合成filterBlock返回 */
  std::string result_;           // Filter data computed so far
  std::vector<Slice> tmp_keys_;  // policy_->CreateFilter() argument
  /** 各个filter在result_中的offset(Starting index in result_ of each filter) */
  std::vector<uint32_t> filter_offsets_;
};

class FilterBlockReader {
 public:
  // REQUIRES: "contents" and *policy must stay live while *this is live.
  FilterBlockReader(const FilterPolicy* policy, const Slice& contents);
  bool KeyMayMatch(uint64_t block_offset, const Slice& key);

 private:
  const FilterPolicy* policy_;
  const char* data_;    // Pointer to filter data (at block-start)
  const char* offset_;  // Pointer to beginning of offset array (at block-end)
  size_t num_;          // Number of entries in offset array
  size_t base_lg_;      // Encoding parameter (see kFilterBaseLg in .cc file)
};

}  // namespace leveldb

#endif  // STORAGE_LEVELDB_TABLE_FILTER_BLOCK_H_
