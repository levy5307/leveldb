// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#ifndef STORAGE_LEVELDB_TABLE_BLOCK_H_
#define STORAGE_LEVELDB_TABLE_BLOCK_H_

#include <stddef.h>
#include <stdint.h>

#include "leveldb/iterator.h"

namespace leveldb {

struct BlockContents;
class Comparator;

/**
 * sstable中的数据以Block为单位存储
 * Block格式：
 *  ________________________________________________________
 * |    block data    |       type      |       crc32       |
 *  --------------------------------------------------------
 * Block data格式：
 *  ____________________________________________________________________________
 * | c1 c2 | c3 c4 | c5 c6 | c1 restart | c3 restart | c5 restart | restart num |
 *  ----------------------------------------------------------------------------
 * c(n) restart代表c(n)的偏移。需要使用偏移是因为c1/c2/c3等长度不一
 * restart num代表restart的数量或者crc信息。上图记录了restart数量
 * block_restart_interval指定了一个分区里有多少个entry，比如上图中就是2
 *
 * c(n) entry格式：
 *  _________________________________________________________________
 * | key共享长度 | key非共享长度 | value长度 | key非共享内容 | value内容 |
 *  -----------------------------------------------------------------
 **/
class Block {
 public:
  // Initialize the block with the specified contents.
  explicit Block(const BlockContents& contents);

  Block(const Block&) = delete;
  Block& operator=(const Block&) = delete;

  ~Block();

  size_t size() const { return size_; }
  Iterator* NewIterator(const Comparator* comparator);

 private:
  class Iter;

  uint32_t NumRestarts() const;

  const char* data_;
  size_t size_;
  uint32_t restart_offset_;  // Offset in data_ of restart array
  /** data_是否自己拥有，也就是说data_字段是否需要自己释放 */
  bool owned_;               // Block owns data_[]
};

}  // namespace leveldb

#endif  // STORAGE_LEVELDB_TABLE_BLOCK_H_
