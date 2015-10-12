#pragma once
#include <iostream>
#include <stdexcept>
#include <list>
#include <string>

class TableFormatStream
{
private:
  std::ostream & output_stream_;
  std::list<std::string> column_names_;
  size_t cursor_;
public:
  template <typename ...Args>
  TableFormatStream(std::ostream & output_stream,
                    const std::string& column_name,
                    Args ... args)
    : TableFormatStream(output_stream, args...)
  {
    column_names_.push_front(column_name);
  }

  TableFormatStream(std::ostream & output_stream)
    : output_stream_(output_stream)
    , column_names_()
    , cursor_(0)
  {
  }

  ~TableFormatStream()
  {
    if (cursor_ != 0) {
      output_stream_ << std::endl;
    }
  }

  void add_column(const std::string & column_name)
  {
    column_names_.push_back(column_name);
  }

  void print_header(const std::string & prefix = "# ")
  {
    output_stream_ << prefix;
    auto iter = column_names_.cbegin();
    if (iter != column_names_.cend()) { // not empty
      output_stream_ << *iter;
      ++iter;
      for (; iter != column_names_.cend(); ++iter) {
        output_stream_ << "\t" << *iter;
      }
    }
    output_stream_ << std::endl;
  }

  void linebreak()
  {
    output_stream_ << std::endl;
    cursor_ = 0;
  }

  template <typename T>
  TableFormatStream & operator<<(const T& obj)
  {
    size_t n_column = column_names_.size();
    if (cursor_ == 0) {
      output_stream_ << obj;
      ++cursor_;
    } else if (cursor_ < n_column) {
      output_stream_ << "\t" << obj;
      ++cursor_;
    } else {
      throw std::domain_error("TableFormatStream::operator<<() : Last character should be std::endl");
    }
    return *this;
  }
};

