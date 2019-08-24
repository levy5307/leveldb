// Copyright (c) 2012 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#include "table/filter_block.h"

#include "leveldb/filter_policy.h"
#include "util/coding.h"

namespace leveldb {

// See doc/table_format.md for an explanation of the filter block format.

// Generate new filter every 2KB of data
static const size_t kFilterBaseLg = 11;
static const size_t kFilterBase = 1 << kFilterBaseLg;

FilterBlockBuilder::FilterBlockBuilder(const FilterPolicy* policy)
    : policy_(policy) {}

void FilterBlockBuilder::StartBlock(uint64_t block_offset) {
  uint64_t filter_index = (block_offset / kFilterBase);
  assert(filter_index >= filter_offsets_.size());
  while (filter_index > filter_offsets_.size()) {
    GenerateFilter();
  }
}

/** 加入key */
void FilterBlockBuilder::AddKey(const Slice& key) {
  Slice k = key;
  start_.push_back(keys_.size());
  keys_.append(k.data(), k.size());
}

/**
 * 把数据拼装成meta block并返回:
 * result最终格式(也就是meta block的格式)：
 *  |      filter data 1
 *  |      filter data 2
 *  |           ...
 *  |      filter data n
 *  |     filter offset 1
 *  |     filter offset 2
 *  |           ...
 *  |     filter offset n
 *  |  beginning of filter offset
 *  V          base
 **/
Slice FilterBlockBuilder::Finish() {
  /** 先根据keys_中的所有key生成一个filter block */
  if (!start_.empty()) {
    GenerateFilter();
  }

  // Append array of per-filter offsets
  /**
   * 在执行Finish之前，result_中就已经存储了各个生成的filter data
   * append各个filter offset到result_中
   **/
  const uint32_t array_offset = result_.size();
  for (size_t i = 0; i < filter_offsets_.size(); i++) {
    PutFixed32(&result_, filter_offsets_[i]);
  }

  /** Append result_加入所有filter offset之前的长度（即：所有filter block的长度和）*/
  PutFixed32(&result_, array_offset);

  /** Append base */
  result_.push_back(kFilterBaseLg);  // Save encoding parameter in result
  return Slice(result_);
}

/**
 * 根据当前的所有keys_，生成一个filter data，追加到result_中，并更新filter_offsets_
 **/
void FilterBlockBuilder::GenerateFilter() {
  const size_t num_keys = start_.size();
  if (num_keys == 0) {
    // Fast path if there are no keys for this filter
    /** 快速方式，也可以删掉这里的分支判断，直接走下面的流程，逻辑相同 */
    filter_offsets_.push_back(result_.size());
    return;
  }

  // Make list of keys from flattened key structure
  start_.push_back(keys_.size());  // Simplify length computation
  /** temp_keys是从keys_中，根据start_（即各个key的偏移）计算得来的keys数组 */
  tmp_keys_.resize(num_keys);
  for (size_t i = 0; i < num_keys; i++) {
    const char* base = keys_.data() + start_[i];
    size_t length = start_[i + 1] - start_[i];
    tmp_keys_[i] = Slice(base, length);
  }

  // Generate filter for current set of keys and append to result_.
  /** 获取当前filter data的offset，存入filter_offsets_中 */
  filter_offsets_.push_back(result_.size());

  /** 根据num_keys个key创建一个filter data，追加到result_中 */
  policy_->CreateFilter(&tmp_keys_[0], static_cast<int>(num_keys), &result_);

  tmp_keys_.clear();
  keys_.clear();
  start_.clear();
}

FilterBlockReader::FilterBlockReader(const FilterPolicy* policy,
                                     const Slice& contents)
    : policy_(policy), data_(nullptr), offset_(nullptr), num_(0), base_lg_(0) {
  size_t n = contents.size();
  if (n < 5) return;  // 1 byte for base_lg_ and 4 for start of offset array
  /** base_lg_ = 途中的base处保存的值 */
  base_lg_ = contents[n - 1];
  /** 这里的last_word = beginning of filter offset处保存的值，data_ + last_word = filter offset 1的地址 */
  uint32_t last_word = DecodeFixed32(contents.data() + n - 5);
  if (last_word > n - 5) return;
  data_ = contents.data();
  offset_ = data_ + last_word;
  /** 通过计算filter offset的数量, 得到filter的数量 */
  num_ = (n - 5 - last_word) / 4;
}

bool FilterBlockReader::KeyMayMatch(uint64_t block_offset, const Slice& key) {
  /** 先根据block_offset和base_lg_解析出filter block的index */
  uint64_t index = block_offset >> base_lg_;
  if (index < num_) {
    /** 解析得到filter block的开始offset和结束offset */
    uint32_t start = DecodeFixed32(offset_ + index * 4);
    uint32_t limit = DecodeFixed32(offset_ + index * 4 + 4);
    /** limit必须小于filter offset(因为不管是start还是limit都是在filter block区间的，不可能超越filter offset这个区间) */
    if (start <= limit && limit <= static_cast<size_t>(offset_ - data_)) {
      Slice filter = Slice(data_ + start, limit - start);
      return policy_->KeyMayMatch(key, filter);
    } else if (start == limit) {
      // Empty filters do not match any keys
      return false;
    }
  }
  return true;  // Errors are treated as potential matches
}

}  // namespace leveldb
