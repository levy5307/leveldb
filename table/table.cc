// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#include "leveldb/table.h"

#include "leveldb/cache.h"
#include "leveldb/comparator.h"
#include "leveldb/env.h"
#include "leveldb/filter_policy.h"
#include "leveldb/options.h"
#include "table/block.h"
#include "table/filter_block.h"
#include "table/format.h"
#include "table/two_level_iterator.h"
#include "util/coding.h"

namespace leveldb {
/**
 * SSTable文件格式：
 *        Data Block 1
 *        Data Block 2        数据
 *          ......      __________________________
 *        Meta Block
 *     Meta Index Block
 *       Index Block
 *          Footer
 *
 *  Data Block/Meta Index Block/Index Block三者都是用block来存储的，所以是利用(key, value)的格式来存储一条条的record的;
 *  而Meta Block则不同
 *
 *  1.Data Block中的KV记录是按照从小到大排序的
 *    Data Block格式：
 *      ________________________________________________________
 *     |    Block Data    |       type      |       crc32       |
 *      --------------------------------------------------------
 *    Block data格式：
 *      ____________________________________________________________________________
 *      | c1 c2 | c3 c4 | c5 c6 | c1 restart | c3 restart | c5 restart | restart num |
 *      ----------------------------------------------------------------------------
 *    c(n) restart代表c(n)的偏移。需要使用偏移是因为c1/c2/c3等长度不一
 *    restart num代表restart的数量或者crc信息。上图记录了restart数量
 *    block_restart_interval指定了一个分区里有多少个entry，比如上图中就是2
 *
 *    c(n) entry格式：
 *      _________________________________________________________________
 *     | key共享长度 | key非共享长度 | value长度 | key非共享内容 | value内容 |
 *      -----------------------------------------------------------------
 *
 *  2.Meta Block就是filter block, 存储了block的filter数据, 用于加快查询的速度
 *    Meta Block格式:
 *     |      filter data 1
 *     |      filter data 2
 *     |           ...
 *     |      filter data n
 *     |     filter offset 1
 *     |     filter offset 2
 *     |           ...
 *     |     filter offset n
 *     |  beginning of filter offset
 *     V          base
 *
 *  3.Meta Index Block:
 *      key: name of meta block i
 *      value: block handle of meta block i
 *
 *  4.Index Block内的每条记录记录是对某个Data Block建立的索引信息，每条索引信息包括三个内容：
 *      key: last key of data block i <= key < first key of data block i + 1
 *      value: block handle of data block i
 *
 *  5.Footer的格式：
 *      offset    size        (metaindex_handle: meta index block handle)
 *      offset    size        (index_handle: index block handle)
 *          padding
 *          magic             (8Bytes litten-endian)
 *
 **/

struct Table::Rep {
  ~Rep() {
    delete filter;
    delete[] filter_data;
    delete index_block;
  }

  Options options;
  Status status;
  RandomAccessFile* file;
  uint64_t cache_id;
  FilterBlockReader* filter;
  const char* filter_data;

  BlockHandle metaindex_handle;  // Handle to metaindex_block: saved from footer
  Block* index_block;
};

Status Table::Open(const Options& options, RandomAccessFile* file,
                   uint64_t size, Table** table) {
  *table = nullptr;
  if (size < Footer::kEncodedLength) {
    return Status::Corruption("file is too short to be an sstable");
  }

  /** 读取footer部分内容 */
  char footer_space[Footer::kEncodedLength];
  Slice footer_input;
  Status s = file->Read(size - Footer::kEncodedLength, Footer::kEncodedLength,
                        &footer_input, footer_space);
  if (!s.ok()) return s;

  /** 读取footer部分内容到Footer里 */
  Footer footer;
  s = footer.DecodeFrom(&footer_input);
  if (!s.ok()) return s;

  // Read the index block
  /** 读取index block的内容 */
  BlockContents index_block_contents;
  if (s.ok()) {
    ReadOptions opt;
    if (options.paranoid_checks) {
      opt.verify_checksums = true;
    }
    s = ReadBlock(file, opt, footer.index_handle(), &index_block_contents);
  }

  if (s.ok()) {
    // We've successfully read the footer and the index block: we're
    // ready to serve requests.
    /** 创建index block */
    Block* index_block = new Block(index_block_contents);
    Rep* rep = new Table::Rep;
    rep->options = options;
    rep->file = file;
    rep->metaindex_handle = footer.metaindex_handle();
    rep->index_block = index_block;
    rep->cache_id = (options.block_cache ? options.block_cache->NewId() : 0);
    rep->filter_data = nullptr;
    rep->filter = nullptr;

    /** 创建table并读取meta */
    *table = new Table(rep);
    (*table)->ReadMeta(footer);
  }

  return s;
}

void Table::ReadMeta(const Footer& footer) {
  if (rep_->options.filter_policy == nullptr) {
    return;  // Do not need any metadata
  }

  // TODO(sanjay): Skip this if footer.metaindex_handle() size indicates
  // it is an empty block.
  ReadOptions opt;
  if (rep_->options.paranoid_checks) {
    opt.verify_checksums = true;
  }

  /** 获取meta index block */
  BlockContents contents;
  if (!ReadBlock(rep_->file, opt, footer.metaindex_handle(), &contents).ok()) {
    // Do not propagate errors since meta info is not needed for operation
    return;
  }
  Block* meta = new Block(contents);

  /** 从meta index block中定位到key, 其value是meta block handle */
  Iterator* iter = meta->NewIterator(BytewiseComparator());
  std::string key = "filter.";
  key.append(rep_->options.filter_policy->Name());
  iter->Seek(key);

  /** 读取该filter block(meta block) */
  if (iter->Valid() && iter->key() == Slice(key)) {
    ReadFilter(iter->value());
  }
  delete iter;
  delete meta;
}

void Table::ReadFilter(const Slice& filter_handle_value) {
  /** 解析获取filter handle(meta block handle) */
  Slice v = filter_handle_value;
  BlockHandle filter_handle;
  if (!filter_handle.DecodeFrom(&v).ok()) {
    return;
  }

  // We might want to unify with ReadBlock() if we start
  // requiring checksum verification in Table::Open.
  ReadOptions opt;
  if (rep_->options.paranoid_checks) {
    opt.verify_checksums = true;
  }

  /** 通过filter handle读取到meta block */
  BlockContents block;
  if (!ReadBlock(rep_->file, opt, filter_handle, &block).ok()) {
    return;
  }
  if (block.heap_allocated) {
    rep_->filter_data = block.data.data();  // Will need to delete later
  }

  /** 返回filter(meta) block reader */
  rep_->filter = new FilterBlockReader(rep_->options.filter_policy, block.data);
}

Table::~Table() { delete rep_; }

static void DeleteBlock(void* arg, void* ignored) {
  delete reinterpret_cast<Block*>(arg);
}

static void DeleteCachedBlock(const Slice& key, void* value) {
  Block* block = reinterpret_cast<Block*>(value);
  delete block;
}

static void ReleaseBlock(void* arg, void* h) {
  Cache* cache = reinterpret_cast<Cache*>(arg);
  Cache::Handle* handle = reinterpret_cast<Cache::Handle*>(h);
  cache->Release(handle);
}

// Convert an index iterator value (i.e., an encoded BlockHandle)
// into an iterator over the contents of the corresponding block.
/** 根据index_value从table中找到对应的block，然后返回该block的iterator(该iterator就是block reader) */
Iterator* Table::BlockReader(void* arg, const ReadOptions& options,
                             const Slice& index_value) {
  Table* table = reinterpret_cast<Table*>(arg);
  Cache* block_cache = table->rep_->options.block_cache;
  Block* block = nullptr;
  Cache::Handle* cache_handle = nullptr;

  BlockHandle handle;
  Slice input = index_value;
  Status s = handle.DecodeFrom(&input);
  // We intentionally allow extra stuff in index_value so that we
  // can add more features in the future.

  if (s.ok()) {
    BlockContents contents;
    /** block_cache不为空的话，从block_cache中获取到对应的block，如果获取不到，则从文件中读 */
    if (block_cache != nullptr) {
      char cache_key_buffer[16];
      EncodeFixed64(cache_key_buffer, table->rep_->cache_id);
      EncodeFixed64(cache_key_buffer + 8, handle.offset());
      Slice key(cache_key_buffer, sizeof(cache_key_buffer));
      cache_handle = block_cache->Lookup(key);
      if (cache_handle != nullptr) {
        /** 从cache中找到了 */
        block = reinterpret_cast<Block*>(block_cache->Value(cache_handle));
      } else {
        /** cache中没有找到，则读取文件，并将读取出来的block加入到block_cache中 */
        s = ReadBlock(table->rep_->file, options, handle, &contents);
        if (s.ok()) {
          block = new Block(contents);
          if (contents.cachable && options.fill_cache) {
            cache_handle = block_cache->Insert(key, block, block->size(),
                                               &DeleteCachedBlock);
          }
        }
      }
    } else {
      /** block_cache为空，则直接从文件中读 */
      s = ReadBlock(table->rep_->file, options, handle, &contents);
      if (s.ok()) {
        block = new Block(contents);
      }
    }
  }

  /** 返回block的iterator */
  Iterator* iter;
  if (block != nullptr) {
    iter = block->NewIterator(table->rep_->options.comparator);
    if (cache_handle == nullptr) {
      iter->RegisterCleanup(&DeleteBlock, block, nullptr);
    } else {
      iter->RegisterCleanup(&ReleaseBlock, block_cache, cache_handle);
    }
  } else {
    iter = NewErrorIterator(s);
  }
  return iter;
}

/** 返回Table的two level iterator */
Iterator* Table::NewIterator(const ReadOptions& options) const {
  return NewTwoLevelIterator(
      rep_->index_block->NewIterator(rep_->options.comparator),
      &Table::BlockReader, const_cast<Table*>(this), options);
}

Status Table::InternalGet(const ReadOptions& options, const Slice& k, void* arg,
                          void (*handle_result)(void*, const Slice&,
                                                const Slice&)) {
  Status s;

  /** 从index block中找到对应的meta block, meta block中存储了对应data block的用于快速查找的filter */
  Iterator* iiter = rep_->index_block->NewIterator(rep_->options.comparator);
  iiter->Seek(k);   // 根据key找到对应的pair(key: key, value: block handle of data clock)
  if (iiter->Valid()) {
    /** block handle of data block */
    Slice handle_value = iiter->value();

    /** 获取到meta block */
    FilterBlockReader* filter = rep_->filter;
    BlockHandle handle;

    /** 先从meta block中查找，如果没有，则一定没有；如果能找到，则不一定有，需要去block再查找验证 */
    if (filter != nullptr && handle.DecodeFrom(&handle_value).ok() &&
        /** 根据data block offset可以获取到meta block index */
        !filter->KeyMayMatch(handle.offset(), k)) {
        /** 若meta block中没有，一定没有 */
      // Not found
    } else {
      /** 若meta block中没有，则不一定有(bloom filter的特性决定的)，则需要从block中重新查找来最终确定是否存在 */
      Iterator* block_iter = BlockReader(this, options, iiter->value());
      block_iter->Seek(k);

      /** 寻找到了，调用传入的回调函数 */
      if (block_iter->Valid()) {
        (*handle_result)(arg, block_iter->key(), block_iter->value());
      }
      s = block_iter->status();
      delete block_iter;
    }
  }
  if (s.ok()) {
    s = iiter->status();
  }
  delete iiter;
  return s;
}

uint64_t Table::ApproximateOffsetOf(const Slice& key) const {
  Iterator* index_iter =
      rep_->index_block->NewIterator(rep_->options.comparator);
  index_iter->Seek(key);
  uint64_t result;
  if (index_iter->Valid()) {
    BlockHandle handle;
    Slice input = index_iter->value();
    Status s = handle.DecodeFrom(&input);
    if (s.ok()) {
      result = handle.offset();
    } else {
      // Strange: we can't decode the block handle in the index block.
      // We'll just return the offset of the metaindex block, which is
      // close to the whole file size for this case.
      result = rep_->metaindex_handle.offset();
    }
  } else {
    // key is past the last key in the file.  Approximate the offset
    // by returning the offset of the metaindex block (which is
    // right near the end of the file).
    result = rep_->metaindex_handle.offset();
  }
  delete index_iter;
  return result;
}

}  // namespace leveldb
