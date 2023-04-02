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

#include <list>
#include <memory>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

namespace IO {

class JsonHandler {
 private:
  nlohmann::json data;

 public:

  /**
   * @brief Construct a new JsonHandler object
   *
   * @param json_encoded_string JSON encoded string
   */
  JsonHandler(std::string_view json_encoded_string);

  /**
   * @brief Construct a new JsonHandler object
   *
   * @param data simdjson data
   */
  JsonHandler(nlohmann::json &data);

  /**
   * @brief Destroy the JsonHandler object
   *
   */
  ~JsonHandler();

  /**
   * @brief Returns the lengths of the data if it can be cast to an array
   *
   * @return int Length of the data
   */
  int length();

  /**
   * @brief Check if JSON data has a certain key
   *
   * @param key
   * @return true If key exists
   * @return false If key doesn't exist
   */
  bool has_key(std::string_view key);

  /**
   * @brief Get the inner JSON object at a string key
   *
   * @param key JSON key
   * @return JsonHandler Inner JSON object at the key
   */
  JsonHandler operator[](std::string_view key);

  /**
   * @brief Get the inner JSON object at an integer key
   *
   * @param index JSON index
   * @return JsonHandler
   */
  JsonHandler operator[](int index);

  /**
   * @brief Get boolean parameter
   *
   * @param key Key of the parameter
   * @return bool Value of the parameter
   */
  bool get_bool(std::string_view key);

  /**
   * @brief Get boolean parameter with default value
   *
   * @param key Key of the parameter
   * @param default_value Value to return if key doesn't exist
   * @return bool Value of the parameter or default value if key doesn't exist
   */
  bool get_bool(std::string_view key, bool default_value);

  /**
   * @brief Get integer parameter
   *
   * @param key Key of the parameter
   * @return int Value of the parameter
   */
  int get_int(std::string_view key);

  /**
   * @brief Get integer parameter with default value
   *
   * @param key Key of the parameter
   * @param default_value Value to return if key doesn't exist
   * @return int Value of the parameter or default value if key doesn't exist
   */
  int get_int(std::string_view key, int default_value);

  /**
   * @brief Get double parameter
   *
   * @param key Key of the parameter
   * @return double Value of the parameter
   */
  double get_double(std::string_view key);

  /**
   * @brief Get double parameter with default value
   *
   * @param key Key of the parameter
   * @param default_value Value to return if key doesn't exist
   * @return double Value of the parameter or default value if key doesn't exist
   */
  double get_double(std::string_view key, double default_value);

  /**
   * @brief Get double array parameter
   *
   * @param key Key of the parameter
   * @return std::vector<double> Value of the parameter
   */
  std::vector<double> get_double_array(std::string_view key);

  /**
   * @brief Get double array parameter with default value
   *
   * @param key Key of the parameter
   * @param default_value Value to return if key doesn't exist
   * @return std::vector<double> Value of the parameter or default value if key
   * doesn't exist
   */
  std::vector<double> get_double_array(std::string_view key,
                                       std::vector<double> default_value);

  /**
   * @brief Get integer array parameter
   *
   * @param key Key of the parameter
   * @return std::vector<int> Value of the parameter
   */
  std::vector<int> get_int_array(std::string_view key);

  /**
   * @brief Get string parameter
   *
   * @param key Key of the parameter
   * @return std::string_view Value of the parameter
   */
  std::string get_string(std::string_view key);

  /**
   * @brief Get string array parameter
   *
   * @param key Key of the parameter
   * @return std::vector<std::string_view> Value of the parameter
   */
  std::vector<std::string> get_string_array(std::string_view key);
};

JsonHandler::JsonHandler(std::string_view json_encoded_string) {
  data = nlohmann::json::parse(json_encoded_string);
}

JsonHandler::JsonHandler(nlohmann::json& data) {
  this->data = data;
}

JsonHandler::~JsonHandler() {}

int JsonHandler::length() {
  return data.size();
}

bool JsonHandler::has_key(std::string_view key) {
  return data.contains(key);
}

JsonHandler JsonHandler::operator[](std::string_view key) {
  return JsonHandler(data[key]);
}

JsonHandler JsonHandler::operator[](int index) {
  return JsonHandler(data[index]);
}

bool JsonHandler::get_bool(std::string_view key) {
  return data[key];
}

bool JsonHandler::get_bool(std::string_view key, bool default_value) {
  return data.value(key, default_value);
}

int JsonHandler::get_int(std::string_view key) {
  return data[key];
}

int JsonHandler::get_int(std::string_view key, int default_value) {
  return data.value(key, default_value);
}

double JsonHandler::get_double(std::string_view key) {
  return data[key];
}

double JsonHandler::get_double(std::string_view key, double default_value) {
  return data.value(key, default_value);
}

std::vector<double> JsonHandler::get_double_array(std::string_view key) {
  std::vector<double> vector;
  if (!data[key].is_array())
  {
    return {data[key]};
  }
  return data[key].get<std::vector<double>>();
}

std::vector<double> JsonHandler::get_double_array(
    std::string_view key, std::vector<double> default_value) {
  if (!has_key(key))
  {
    return default_value;
  }
  if (!data[key].is_array())
  {
    return {data[key]};
  }
  return data[key].get<std::vector<double>>();
}

std::vector<int> JsonHandler::get_int_array(std::string_view key) {
  if (!data[key].is_array())
  {
    return {data[key]};
  }
  return data[key].get<std::vector<int>>();
}

std::string JsonHandler::get_string(std::string_view key) {
  return data[key];
}

std::vector<std::string> JsonHandler::get_string_array(
    std::string_view key) {
  return data[key].get<std::vector<std::string>>();
}

}  // namespace IO

#endif  // SVZERODSOLVER_IO_JSONHANDLER_HPP_
