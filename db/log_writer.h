// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#ifndef STORAGE_LEVELDB_DB_LOG_WRITER_H_
#define STORAGE_LEVELDB_DB_LOG_WRITER_H_

#include <stdint.h>

#include "db/log_format.h"
#include "leveldb/slice.h"
#include "leveldb/status.h"

namespace leveldb {

class WritableFile;
/**
 * log文件中的数据是以block为单位组织的。写日志时，基于一致性考虑，并没有按block单位写。
 * 每次更新均对log文件进行IO，根据WritOption::sync决定是否做强制sync，读取时以block为单位做IO以及校验
 * log文件格式：
 *  init data
 *  block 1
 *  block 2
 *   ....
 *  block n
 * 当前版本中init data为空，block中存储实际的数据
 * block的整体结构如下：
 *      record 0
 *      record 1
 *       ....
 *      record n
 *      tailer
 * tailer: 如果block最后剩余部分小于record的头部，则将其填充为0作为tailer。而record写入下一个block
 * record: 每次写入作为一个record, record组成：
 *     (地址)低--->高：
 *      ________________________________________________________________
 *     | checksum(uint32) | length(uint16) | type(uint8) | data(length) |
 *      ----------------------------------------------------------------
 *     checksum: type和data的crc校验
 *     length: data的长度
 *     type: 为了避免block内部碎片的产生，一份record可能会跨block，所以根据record内保存数据占更新写入数据的完整与否，
 *           当前分为4种type：FULL, FIRST, MIDDLE, LAST, 依次表示record内保存的是完整数据的全部、开始、中间或最后部分
 *
 * log的写入时顺序写，读取只会在启动时发生，不会是性能瓶颈，log的数据也就没有进行压缩
 **/

namespace log {

class Writer {
 public:
  // Create a writer that will append data to "*dest".
  // "*dest" must be initially empty.
  // "*dest" must remain live while this Writer is in use.
  explicit Writer(WritableFile* dest);

  // Create a writer that will append data to "*dest".
  // "*dest" must have initial length "dest_length".
  // "*dest" must remain live while this Writer is in use.
  Writer(WritableFile* dest, uint64_t dest_length);

  Writer(const Writer&) = delete;
  Writer& operator=(const Writer&) = delete;

  ~Writer();

  Status AddRecord(const Slice& slice);

 private:
  Status EmitPhysicalRecord(RecordType type, const char* ptr, size_t length);

  WritableFile* dest_;
  int block_offset_;  // Current offset in block

  // crc32c values for all supported record types.  These are
  // pre-computed to reduce the overhead of computing the crc of the
  // record type stored in the header.
  uint32_t type_crc_[kMaxRecordType + 1];
};

}  // namespace log
}  // namespace leveldb

#endif  // STORAGE_LEVELDB_DB_LOG_WRITER_H_
