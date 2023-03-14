// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/**
 * @file jsonhandler.hpp
 * @brief IO::JsonHandler source file
 */
#ifndef SVZERODSOLVER_IO_JSONHANDLER_HPP_
#define SVZERODSOLVER_IO_JSONHANDLER_HPP_

#include <stdexcept>
#include <string>

#include "simdjson.h"

namespace IO {

class JsonHandler {
 private:
  simdjson::simdjson_result<simdjson::dom::element> data;
  simdjson::dom::parser *parser;

 public:
  JsonHandler(std::string_view json_encoded_string);
  JsonHandler(simdjson::simdjson_result<simdjson::dom::element> data);
  ~JsonHandler();

  int length();
  bool has_key(std::string_view key);

  JsonHandler operator[](std::string_view key);
  JsonHandler operator[](int index);

  bool get_bool(std::string_view key);
  bool get_bool(std::string_view key, bool default_value);

  int get_int(std::string_view key);
  int get_int(std::string_view key, int default_value);

  double get_double(std::string_view key);
  double get_double(std::string_view key, double default_value);

  std::vector<double> get_double_array(std::string_view key);
  std::vector<double> get_double_array(std::string_view key,
                                       std::vector<double> default_value);

  std::vector<int> get_int_array(std::string_view key);

  std::string_view get_string(std::string_view key);

  std::vector<std::string_view> get_string_array(std::string_view key);

  std::list<JsonHandler> get_list(std::string_view key);
};

JsonHandler::JsonHandler(std::string_view json_encoded_string) {
  parser = new simdjson::dom::parser();
  auto string = simdjson::padded_string(json_encoded_string);
  data = parser->parse(string);
}
JsonHandler::JsonHandler(
    simdjson::simdjson_result<simdjson::dom::element> data) {
  this->data = data;
}

JsonHandler::~JsonHandler() {}

int JsonHandler::length() {
  int size;
  try {
    size = data.get_array().size();
  } catch (simdjson::simdjson_error) {
    throw std::runtime_error("Not an array");
  }
  return size;
}

bool JsonHandler::has_key(std::string_view key) {
  bool has_key = false;
  for (auto &&field : data.get_object()) {
    if (field.key == key) {
      has_key = true;
    }
  }
  return has_key;
}

JsonHandler JsonHandler::operator[](std::string_view key) {
  if (!has_key(key)) {
    throw std::runtime_error("Key not found: " + static_cast<std::string>(key));
  }
  return JsonHandler(data.at_key(key));
}

JsonHandler JsonHandler::operator[](int index) {
  simdjson::simdjson_result<simdjson::dom::element> output;
  try {
    output = data.get_array().at(index);
  } catch (simdjson::simdjson_error) {
    throw std::runtime_error("Index out of range: " + std::to_string(index));
  }
  return JsonHandler(output);
}

bool JsonHandler::get_bool(std::string_view key) {
  bool output;
  try {
    output = data.at_key(key).get_bool();
  } catch (simdjson::simdjson_error) {
    throw std::runtime_error("Key not found: " + static_cast<std::string>(key));
  }
  return output;
}

bool JsonHandler::get_bool(std::string_view key, bool default_value) {
  try {
    return data.at_key(key).get_bool();
  } catch (simdjson::simdjson_error) {
    return default_value;
  }
}

int JsonHandler::get_int(std::string_view key) {
  int output;
  try {
    output = data.at_key(key).get_int64();
  } catch (simdjson::simdjson_error) {
    throw std::runtime_error("Key not found: " + static_cast<std::string>(key));
  }
  return output;
}

int JsonHandler::get_int(std::string_view key, int default_value) {
  try {
    return data.at_key(key).get_int64();
  } catch (simdjson::simdjson_error) {
    return default_value;
  }
}

double JsonHandler::get_double(std::string_view key) {
  double output;
  try {
    output = data.at_key(key).get_double();
  } catch (simdjson::simdjson_error) {
    throw std::runtime_error("Key not found: " + static_cast<std::string>(key));
  }
  return output;
}

double JsonHandler::get_double(std::string_view key, double default_value) {
  try {
    return data.at_key(key).get_double();
  } catch (simdjson::simdjson_error) {
    return default_value;
  }
}

std::vector<double> JsonHandler::get_double_array(std::string_view key) {
  std::vector<double> vector;
  try {
    for (auto x : data.at_key(key).get_array()) {
      vector.push_back(x.get_double());
    }
  } catch (simdjson::simdjson_error) {
    try {
      vector.push_back(data.at_key(key).get_double());
    } catch (simdjson::simdjson_error) {
      throw std::runtime_error("Key not found: " +
                               static_cast<std::string>(key));
    }
  };
  return vector;
}

std::vector<double> JsonHandler::get_double_array(
    std::string_view key, std::vector<double> default_value) {
  std::vector<double> vector;
  try {
    for (auto x : data.at_key(key).get_array()) {
      vector.push_back(x.get_double());
    }
  } catch (simdjson::simdjson_error) {
    try {
      vector.push_back(data.at_key(key).get_double());
    } catch (simdjson::simdjson_error) {
      return default_value;
    }
  };
  return vector;
}

std::vector<int> JsonHandler::get_int_array(std::string_view key) {
  std::vector<int> vector;
  try {
    for (auto x : data.at_key(key).get_array()) {
      int val = x.get_int64();
      vector.push_back(val);
    }
  } catch (simdjson::simdjson_error) {
    try {
      int val = data.at_key(key).get_int64();
      vector.push_back(val);
    } catch (simdjson::simdjson_error) {
      throw std::runtime_error("Key not found: " +
                               static_cast<std::string>(key));
    }
  };
  return vector;
}

std::string_view JsonHandler::get_string(std::string_view key) {
  std::string_view output;
  try {
    output = data.at_key(key).get_string();
  } catch (simdjson::simdjson_error) {
    throw std::runtime_error("Key not found: " + static_cast<std::string>(key));
  }
  return output;
}

std::vector<std::string_view> JsonHandler::get_string_array(
    std::string_view key) {
  std::vector<std::string_view> vector;
  try {
    for (auto x : data.at_key(key).get_array()) {
      vector.push_back(x.get_string());
    }
  } catch (simdjson::simdjson_error) {
    try {
      vector.push_back(data.at_key(key).get_string());
    } catch (simdjson::simdjson_error) {
      throw std::runtime_error("Key not found: " +
                               static_cast<std::string>(key));
    }
  };
  return vector;
}

}  // namespace IO

#endif  // SVZERODSOLVER_IO_JSONHANDLER_HPP_