// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.
//
// Endian-neutral encoding:
// * Fixed-length numbers are encoded with least-significant byte first
// * In addition we support variable length "varint" encoding
// * Strings are encoded prefixed by their length in varint format

#ifndef STORAGE_LEVELDB_UTIL_CODING_H_
#define STORAGE_LEVELDB_UTIL_CODING_H_

#include <cstdint>
#include <cstring>
#include <string>

#include "leveldb/slice.h"
#include "port/port.h"

namespace leveldb {

/**
 * Fixed: 固定长度数字
 * Varint: 变长数字，每一个byte用低7个字节来表示，第8个字节用于表示是否还有后续的部分, 对于大数，有可能用到5 Bytes
 * */
// Standard Put... routines append to a string
void PutFixed32(std::string* dst, uint32_t value);
void PutFixed64(std::string* dst, uint64_t value);
void PutVarint32(std::string* dst, uint32_t value);
void PutVarint64(std::string* dst, uint64_t value);
void PutLengthPrefixedSlice(std::string* dst, const Slice& value);

// Standard Get... routines parse a value from the beginning of a Slice
// and advance the slice past the parsed value.
bool GetVarint32(Slice* input, uint32_t* value);
bool GetVarint64(Slice* input, uint64_t* value);
bool GetLengthPrefixedSlice(Slice* input, Slice* result);

// Pointer-based variants of GetVarint...  These either store a value
// in *v and return a pointer just past the parsed value, or return
// nullptr on error.  These routines only look at bytes in the range
// [p..limit-1]
const char* GetVarint32Ptr(const char* p, const char* limit, uint32_t* v);
const char* GetVarint64Ptr(const char* p, const char* limit, uint64_t* v);

// Returns the length of the varint32 or varint64 encoding of "v"
int VarintLength(uint64_t v);

// Lower-level versions of Put... that write directly into a character buffer
// and return a pointer just past the last byte written.
// REQUIRES: dst has enough space for the value being written
char* EncodeVarint32(char* dst, uint32_t value);
char* EncodeVarint64(char* dst, uint64_t value);

// TODO(costan): Remove port::kLittleEndian and the fast paths based on
//               std::memcpy when clang learns to optimize the generic code, as
//               described in https://bugs.llvm.org/show_bug.cgi?id=41761
//
// The platform-independent code in DecodeFixed{32,64}() gets optimized to mov
// on x86 and ldr on ARM64, by both clang and gcc. However, only gcc optimizes
// the platform-independent code in EncodeFixed{32,64}() to mov / str.

// Lower-level versions of Put... that write directly into a character buffer
// REQUIRES: dst has enough space for the value being written

/**
 * DecodeFixed中的平台无关代码，对于clang和gcc情况下都做了优化(优化成了单个mov/ldr指令),
 * 但是EncodeFixed中的平台无关代码只有在gcc情况下做了优化(即clang对于store操作的支持不够优秀)
 * 但是两者（DecodeFixed和EncodeFixed中的平台无关代码）对于msvc都没有做优化（msvc是微软的VC编译器）
 * 可以从如下网址查看到编译后形成的汇编代码： https://godbolt.org/z/45S0ID
 **/
inline void EncodeFixed32(char* dst, uint32_t value) {
  uint8_t* const buffer = reinterpret_cast<uint8_t*>(dst);

  if (port::kLittleEndian) {
    // Fast path for little-endian CPUs. All major compilers optimize this to a
    // single mov (x86_64) / str (ARM) instruction.
    std::memcpy(buffer, &value, sizeof(uint32_t));
    return;
  }

  /**
   * 只有gcc情况下对下述平台无关操作做了优化（优化成了单个的mov/ldr指令）, 其他情况下并未做优化
   * 为了尽可能提高效率，才有了上面的对little endian的特殊处理（因为在littlen endian情况下，memcpy做了优化）
   * 否则，就可以直接使用下述平台（big endian/little endian）无关代码
   **/
  // Platform-independent code.
  // Currently, only gcc optimizes this to a single mov / str instruction.
  buffer[0] = static_cast<uint8_t>(value);
  buffer[1] = static_cast<uint8_t>(value >> 8);
  buffer[2] = static_cast<uint8_t>(value >> 16);
  buffer[3] = static_cast<uint8_t>(value >> 24);
}

inline void EncodeFixed64(char* dst, uint64_t value) {
  uint8_t* const buffer = reinterpret_cast<uint8_t*>(dst);

  if (port::kLittleEndian) {
    // Fast path for little-endian CPUs. All major compilers optimize this to a
    // single mov (x86_64) / str (ARM) instruction.
    std::memcpy(buffer, &value, sizeof(uint64_t));
    return;
  }

  // Platform-independent code.
  // Currently, only gcc optimizes this to a single mov / str instruction.
  buffer[0] = static_cast<uint8_t>(value);
  buffer[1] = static_cast<uint8_t>(value >> 8);
  buffer[2] = static_cast<uint8_t>(value >> 16);
  buffer[3] = static_cast<uint8_t>(value >> 24);
  buffer[4] = static_cast<uint8_t>(value >> 32);
  buffer[5] = static_cast<uint8_t>(value >> 40);
  buffer[6] = static_cast<uint8_t>(value >> 48);
  buffer[7] = static_cast<uint8_t>(value >> 56);
}

// Lower-level versions of Get... that read directly from a character buffer
// without any bounds checking.

inline uint32_t DecodeFixed32(const char* ptr) {
  const uint8_t* const buffer = reinterpret_cast<const uint8_t*>(ptr);

  if (port::kLittleEndian) {
    // Fast path for little-endian CPUs. All major compilers optimize this to a
    // single mov (x86_64) / ldr (ARM) instruction.
    uint32_t result;
    std::memcpy(&result, buffer, sizeof(uint32_t));
    return result;
  }

  // Platform-independent code.
  // Clang and gcc optimize this to a single mov / ldr instruction.
  return (static_cast<uint32_t>(buffer[0])) |
         (static_cast<uint32_t>(buffer[1]) << 8) |
         (static_cast<uint32_t>(buffer[2]) << 16) |
         (static_cast<uint32_t>(buffer[3]) << 24);
}

inline uint64_t DecodeFixed64(const char* ptr) {
  const uint8_t* const buffer = reinterpret_cast<const uint8_t*>(ptr);

  if (port::kLittleEndian) {
    // Fast path for little-endian CPUs. All major compilers optimize this to a
    // single mov (x86_64) / ldr (ARM) instruction.
    uint64_t result;
    std::memcpy(&result, buffer, sizeof(uint64_t));
    return result;
  }

  // Platform-independent code.
  // Clang and gcc optimize this to a single mov / ldr instruction.
  return (static_cast<uint64_t>(buffer[0])) |
         (static_cast<uint64_t>(buffer[1]) << 8) |
         (static_cast<uint64_t>(buffer[2]) << 16) |
         (static_cast<uint64_t>(buffer[3]) << 24) |
         (static_cast<uint64_t>(buffer[4]) << 32) |
         (static_cast<uint64_t>(buffer[5]) << 40) |
         (static_cast<uint64_t>(buffer[6]) << 48) |
         (static_cast<uint64_t>(buffer[7]) << 56);
}

// Internal routine for use by fallback path of GetVarint32Ptr
const char* GetVarint32PtrFallback(const char* p, const char* limit,
                                   uint32_t* value);
inline const char* GetVarint32Ptr(const char* p, const char* limit,
                                  uint32_t* value) {
  /** 这里处理的是该数字只有一个字节的情况, 个人感觉这里多此一举，直接调GetVarint32PtrFallback应该就可以了 */
  if (p < limit) {
    uint32_t result = *(reinterpret_cast<const uint8_t*>(p));
    if ((result & 128) == 0) {
      *value = result;
      return p + 1;
    }
  }
  return GetVarint32PtrFallback(p, limit, value);
}

}  // namespace leveldb

#endif  // STORAGE_LEVELDB_UTIL_CODING_H_
