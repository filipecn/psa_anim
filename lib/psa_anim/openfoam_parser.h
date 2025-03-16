/// Copyright (c) 2022, FilipeCN.
///
/// The MIT License (MIT)
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to
/// deal in the Software without restriction, including without limitation the
/// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
/// sell copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in
/// all copies or substantial portions of the Software.
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
/// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
/// IN THE SOFTWARE.
///
///\file openfoam_parser.h
///\author FilipeCN (filipedecn@gmail.com)
///\date 2022-01-11
///
///\brief

#ifndef PSA_ANIM_VIS_OPENFOAM_VIEWER_OPENFOAM_PARSER_H
#define PSA_ANIM_VIS_OPENFOAM_VIEWER_OPENFOAM_PARSER_H

#include <cctype>
#include <hermes/common/file_system.h>

// *********************************************************************************************************************
//                                                                                                       OpenFoamDict
// *********************************************************************************************************************
struct OpenFoamDict {
  // *******************************************************************************************************************
  //                                                                                           OpenFoamDict::DictNode
  // *******************************************************************************************************************
  /// \brief Dictionary Node structure (contains inner dictionary nodes)
  struct DictNode {
    /// Node entry name
    std::string name;
    /// Node value
    std::string value;
    /// Node children
    std::map<std::string, DictNode> fields;
    /// Retrieves a sub-node in the first node level of this dict node
    /// \param field_name node name
    /// \return the dictionary node or an empty node if not found
    const DictNode &operator[](const std::string &field_name) const {
      static DictNode dummy{"not found", {}, {}};
      auto it = fields.find(field_name);
      if (it == fields.end())
        return dummy;
      return it->second;
    }
    /// Dumps the dictionary
    /// \param level tab level for this node
    /// \return dump content
    [[nodiscard]] hermes::Str print(int level = 1) const {
      hermes::Str s;
      std::string base(level * 2, ' ');
      s.appendLine(base, " -> ", value);
      for (const auto &f : fields) {
        s.appendLine(base, f.first);
        s += f.second.print(level + 1);
      }
      return s;
    }
  };
  // *******************************************************************************************************************
  //                                                                                                     CONSTRUCTORS
  // *******************************************************************************************************************
  /// \brief Constructs structure from dictionary path
  /// \param path dictionary path
  explicit OpenFoamDict(const hermes::Path &path) { parseFile(path); }
  // *******************************************************************************************************************
  //                                                                                                        OPERATORS
  // *******************************************************************************************************************
  /// \brief Retrieves a sub-dictionary node
  /// \param field_name node name
  /// \return the dictionary node or an empty node if not found
  const DictNode &operator[](const std::string &field_name) const {
    static DictNode dummy{"none", "not found", {}};
    auto it = dict_.find(field_name);
    if (it == dict_.end())
      return dummy;
    return it->second;
  }
  // *******************************************************************************************************************
  //                                                                                                   STATIC METHODS
  // *******************************************************************************************************************
  /// \brief Parses a string containing a array of values.
  /// Values are expected to be written as n ( v1 v2 v3 ... vn), not necessarily
  /// in the same line.
  /// \tparam T value data type
  /// \param s raw dictionary content
  /// \return a vector of parsed values
  template <typename T>
  static std::vector<T> parseValuesFrom(const std::string &s) {
    if (!(typeid(T) == typeid(double) || typeid(T) == typeid(float))) {
      HERMES_LOG_ERROR("unsupported data type for value list!");
      exit(-1);
    }
    // TODO a field can be uniform, I should check it...
    // TODO for now I'm just assuming the caller of this function will
    // TODO interpret an empty return as a uniform 0
    size_t i = 0;
    std::vector<T> values;
    // read array size
    consumeAll(i, s);
    // if(s == "uniform 0")
    //   return {};
    HERMES_ASSERT(std::isdigit(s[i]));
    auto count_token = nextToken(i, s);
    size_t count = std::stoull(count_token);
    if (count == 0)
      return values;
    // assume values are separated by lines
    while (i < s.size()) {
      consumeAll(i, s);
      auto v = nextToken(i, s);
      HERMES_ASSERT(hermes::Str::isNumber(v));
      T value = 0;
      try {
        value = std::stod(v);
      } catch (std::exception &e) {
        // HERMES_LOG_ERROR(e.what());
        // HERMES_LOG_VARIABLE(v);
        value = 0;
      }
      values.emplace_back(value);
      consumeAll(i, s);
    }
    if (values.size() < count)
      values.resize(count, values.empty() ? 0 : values[0]);
    HERMES_ASSERT(values.size() == count);
    return values;
  }
  /// \brief Parses a string containing a array of lists.
  /// Values are expected to be written as n ( n1(v1 ... vn1) ... nn(v1 ... vnn)
  /// ), not necessarily in the same line.
  /// \tparam T value data type
  /// \param s raw dictionary content
  /// \return a vector of parsed values
  template <typename T>
  static void parseValueListsFrom(const std::string &s,
                                  std::vector<std::vector<T>> &values) {
    HERMES_NOT_IMPLEMENTED
    exit(-1);
    values.clear();
    size_t i = 0;
    while (i < s.size()) {
      auto l = nextLine(i, s);
      auto vss = hermes::Str::split(
          hermes::Str::split(hermes::Str::strip(l, " \n)"), "(")[1]);
      std::vector<T> vs;
      for (const auto &v : vss)
        vs.emplace_back(std::stold(v));
      values.emplace_back(vs);
    }
  }
  // *******************************************************************************************************************
  //                                                                                                          METHODS
  // *******************************************************************************************************************
  /// Parses a file into the dictionary
  /// \param path file path
  void parseFile(const hermes::Path &path) {
    // HERMES_LOG("parsing dict file {}", path);
    HERMES_LOG_AND_RETURN_IF_NOT(
        path.isFile(),
        hermes::Str::concat("file not found! ", path.fullName()).c_str());
    auto s = path.read();
    parse(s);
  }
  ///
  /// \return
  [[nodiscard]] const auto &nodes() const { return dict_; }
  // Prints parsed content
  void print(std::ostream &o = std::cout) {
    o << "OpenFOAM Dictionary [" << dict_.size() << " entries]\n";
    for (const auto &node : dict_)
      o << node.first << ": " << node.second.print();
  }

private:
  /// \brief Moves index i past one block comment (if any).
  /// The first two characters of s must be '/*' or '//'
  /// \note This function jumps only one comment block
  /// \param i starting index
  /// \param s raw string
  static void consumeComment(size_t &i, const std::string &s) {
    // ... things like /*...*//* or /*...*///
    size_t n = s.size();
    if (i < n && s[i] == '/' && i + 1 < n && s[i + 1] == '*')
      for (i++; i < n && !(s[i - 1] == '*' && s[i] == '/'); ++i)
        ;
    else if (i + 1 < n && s[i + 1] == '/')
      for (; i < n && s[i] != '\n'; ++i)
        ;
    else
      return;
    i++;
  }
  /// \brief Moves index i past any sequence of blank characters.
  /// Blank characters are: ' ', '\t', '\n' and ';'
  /// \param i starting index
  /// \param s raw string
  static void consumeBlank(size_t &i, const std::string &s) {
    auto n = s.size();
    size_t j = i;
    do {
      j = i;
      // consume blank
      while (i < n &&
             (s[i] == ' ' || s[i] == '\t' || s[i] == '\n' || s[i] == ';'))
        i++;
    } while (j != i);
  }
  /// \brief Moves index i past any sequence of comments and blank characters.
  /// This function can jump multiple comment blocks and comment lines.
  /// \param i starting index
  /// \param s raw string
  static void consumeAll(size_t &i, const std::string &s) {
    auto n = s.size();
    size_t j = i;
    do {
      j = i;
      consumeComment(i, s);
      consumeBlank(i, s);
    } while (j != i);
  }
  /// \brief Reads a code block node.
  /// \param i receives starting index, and stores end index
  /// \param s raw string
  /// \return string containing the block of code
  static std::string nextCode(size_t &i, const std::string &s) {
    auto n = s.size();
    auto start = i;
    if (s[i] != '#' || i + 1 >= n || s[i + 1] != '{')
      return "";
    ++i;
    while (i < n && !(s[i - 1] == '#' && s[i] == '}'))
      ++i;
    ++i;
    HERMES_ASSERT(s[i] == ';');
    return s.substr(start, i - start);
  }
  /// \brief Read a list of values delimited by '(' and ')'.
  /// \note This function assumes the opening '(' has been previously consumed.
  /// \note The values can be lists.
  /// \param i receives starting index, and stores end index
  /// \param s raw string
  /// \return string containing array
  static std::string nextArray(size_t &i, const std::string &s) {
    consumeAll(i, s);
    auto n = s.size();
    auto start = i;
    int st = 1;
    size_t last_close = 0;
    size_t element_count = 0;
    while (st && i < n) {
      if (s[i] == '(')
        st++;
      if (s[i] == ')') {
        last_close = i;
        st--;
        element_count++;
      }
      i++;
    }
    HERMES_ASSERT(st == 0);
    HERMES_ASSERT(last_close >= start);
    return s.substr(start, last_close - start);
  }
  /// \brief Parses the next token in the string.
  /// \note Consumes any blank character at the beginning of s[i].
  /// \note A token is any sequence of alpha-numeric characters.
  /// \note The token may include special characters '.', '_'
  /// \param i receives starting index, and stores end index
  /// \param s raw string
  /// \return string containing token
  static std::string nextToken(size_t &i, const std::string &s) {
    consumeBlank(i, s);
    auto n = s.size();
    size_t start = i;
    while (i < n && (std::isalnum(s[i]) || s[i] == '.' || s[i] == '_' ||
                     s[i] == '-' || s[i] == '+'))
      ++i;
    if (start == i)
      return "";
    return s.substr(start, (i - 1) - start + 1);
  }
  /// \brief Reads a dict node value.
  /// The value can be a code snippet, a string or a number.
  /// \param i receives starting index, and stores end index
  /// \param s raw string
  /// \return string containing value
  static std::string nextValue(size_t &i, const std::string &s) {
    consumeAll(i, s);
    if (s[i] == '#') // code
      return nextCode(i, s);
    size_t start = i;
    auto n = s.size();
    while (i < n && s[i] != ';')
      ++i;
    auto value = s.substr(start, i - start);
    // before returning, if value contains a list, then we need to process it!
    // TODO improve
    // in order to detect a list, we should expect something like:
    // nonuniform List<*> n ( v1 ... vn )
    size_t j = 0;
    n = value.size();
    auto t = nextToken(j, value);
    if (t == "nonuniform") {
      // if the field is nonuniform, lets really expect a list here!
      consumeAll(j, value);
      // in other words, the next token MUST be List
      t = nextToken(j, value);
      HERMES_ASSERT(t == "List");
      // consume until '>'
      while (j < n && value[j] != '>')
        ++j;
      HERMES_ASSERT(j < n);
      HERMES_ASSERT(value[j] == '>');
      // jump
      ++j;
      consumeAll(j, value);
      // expect simple array now n (v1 ... vn)
      t = nextToken(j, value);
      HERMES_ASSERT(hermes::Str::isInteger(t));
      auto count = t;
      consumeAll(j, value);
      HERMES_ASSERT(value[j] == '(');
      // lets consume it
      j++;
      return count + " " + nextArray(j, value);
    }
    return value;
  }

  static std::string nextLine(size_t &i, const std::string &s) {
    auto n = s.size();
    size_t start = i;
    while (i < n && s[i] != '\n')
      ++i;
    ++i;
    return s.substr(start, (i - 1) - start + 1);
  }
  static std::string nextString(size_t &i, const std::string &s) {
    auto n = s.size();
    size_t start = i;
    size_t end = i;
    // here we need to find a first " and a second "
    bool found_start = false;
    while (i < n) {
      if (s[i] == '\"') {
        if (found_start) {
          end = i++;
          break;
        } else {
          start = i;
          found_start = true;
        }
      }
      ++i;
    }
    HERMES_ASSERT(found_start && start < end);
    return s.substr(start, end - start + 1);
  }
  /// \brief Recursively reads dict nodes delimited by '{' and '}'.
  /// \note This function assumes the opening '{' has been previously consumed.
  /// \param i receives starting index, and stores end index
  /// \param s raw string
  /// \return a set child node structures mapped by their names
  static std::map<std::string, DictNode> readFields(size_t &i,
                                                    const std::string &s) {
    std::map<std::string, DictNode> fields;
    auto n = s.size();
    while (i < n && s[i] != '}') {
      DictNode node;
      auto name = nextToken(i, s);
      i++;
      consumeAll(i, s);
      // the value of a node name can be:
      // - a set of nodes, starting with '{'
      // - a single value
      if (s[i] == '{') {
        i++;
        node.fields = readFields(i, s);
      } else if (s[i] == '\"') {
        // read until next "
        node.value = nextString(i, s);
      } else
        node.value = nextValue(i, s);
      i++;
      consumeAll(i, s);
      node.name = name;
      fields[name] = node;
    }
    i++;
    return fields;
  }
  /// Parses a string into a set of dictionary nodes with names and values
  /// \param s raw text
  void parse(const std::string &s) {
    auto n = s.size();
    size_t i = 0;
    // get rid of any initial comments and blank lines
    consumeAll(i, s);
    while (i < n) {
      // the next character MUST be an letter or number
      HERMES_ASSERT(std::isalnum(s[i]));
      DictNode node;
      auto name = nextToken(i, s);
      // consume any thing between the node name and its value;
      consumeAll(i, s);
      // the value of a node name can be:
      // - a set of nodes, starting with '{'
      // - a list of elements, starting with '('
      // - a single value
      // TODO: should improve this classification, lists might have keywords
      // too! (like List<...)
      if (hermes::Str::isInteger(name)) {
        node.value = name + " ";
        // we have a list!
        consumeAll(i, s);
        HERMES_ASSERT(i < n);
        // the list has the pattern (v1 ... v2)
        // however, v* can be lists too!
        // now name contains the amount of v* in the list.
        // the next character must be '(' or a '{'
        HERMES_ASSERT(s[i] == '(' || s[i] == '{');
        // lets consume it
        i++;
        if (s[i - 1] == '(') {
          node.value += nextArray(i, s);
        } else {
          // TODO here is a little hack, in case of a '{' lets expect a single
          // value in the list!
          consumeAll(i, s);
          node.value += nextToken(i, s);
          consumeAll(i, s);
          HERMES_ASSERT(s[i] == '}');
          ++i;
        }
        name = "list";
      } else if (s[i] == '{') {
        i++;
        node.fields = readFields(i, s);
      } else
        node.value = nextValue(i, s);
      consumeAll(i, s);
      node.name = name;
      dict_[name] = node;
    }
  }

  std::map<std::string, DictNode> dict_;
  //  std::vector<std::string> list;
};
/// \brief Parses a list of 3 component lists of numbers and stores in a vector
/// of vec3. The input string is expected to contain the pattern n (v v v)_1 ...
/// (v v v)_n.
/// \note The input string does not have the main enclosing parenthesis.
/// \note The string can contain spaces, commentaries and break lines.
/// \param s raw string containing the list
/// \return vector of vec3
template <>
inline std::vector<hermes::vec3>
OpenFoamDict::parseValuesFrom(const std::string &s) {
  // HERMES_LOG("Parsing vector of vectors");
  size_t i = 0;
  std::vector<hermes::vec3> values;
  // read array size
  consumeAll(i, s);
  HERMES_ASSERT(std::isdigit(s[i]));
  auto count_token = nextToken(i, s);
  size_t count = std::stoull(count_token);
  if (count == 0)
    return values;
  while (i < s.size()) {
    hermes::vec3 vec;
    consumeAll(i, s);
    // consume opening parenthesis
    HERMES_ASSERT(s[i] == '(');
    ++i;
    for (size_t j = 0; j < 3; ++j) {
      consumeAll(i, s);
      auto v = nextToken(i, s);
      try {
        vec[j] = std::stod(v);
      } catch (std::exception &e) {
        vec[j] = 0;
      }
    }
    consumeAll(i, s);
    // consume closing parenthesis
    HERMES_ASSERT(s[i] == ')');
    ++i;
    values.emplace_back(vec);
    consumeAll(i, s);
  }
  if (values.size() < count)
    values.resize(count, values.empty() ? hermes::vec3{} : values[0]);
  HERMES_ASSERT(values.size() == count);
  return values;
}
/// \brief Parses a list of 3 component lists of numbers and stores in a vector
/// of points. The input string is expected to contain the pattern n (v v v)_1
/// ... (v v v)_n.
/// \note The input string does not have the main enclosing parenthesis.
/// \note The string can contain spaces, commentaries and break lines.
/// \param s raw string containing the list
/// \return vector of points
template <>
inline std::vector<hermes::point3>
OpenFoamDict::parseValuesFrom(const std::string &s) {
  HERMES_LOG("Parsing vector of points");
  size_t i = 0;
  std::vector<hermes::point3> values;
  // read array size
  consumeAll(i, s);
  HERMES_ASSERT(std::isdigit(s[i]));
  auto count_token = nextToken(i, s);
  size_t count = std::stoull(count_token);
  if (count == 0)
    return values;
  while (i < s.size()) {
    hermes::point3 point;
    consumeAll(i, s);
    // consume opening parenthesis
    HERMES_ASSERT(s[i] == '(');
    ++i;
    for (size_t j = 0; j < 3; ++j) {
      consumeAll(i, s);
      auto v = nextToken(i, s);
      try {
        point[j] = std::stof(v);
      } catch (std::exception &e) {
        point[j] = 0;
      }
    }
    consumeAll(i, s);
    // consume closing parenthesis
    HERMES_ASSERT(s[i] == ')');
    ++i;
    values.emplace_back(point);
    consumeAll(i, s);
  }
  if (values.size() < count)
    values.resize(count, values.empty() ? hermes::point3{} : values[0]);
  HERMES_ASSERT(values.size() == count);
  return values;
}
/// \brief Parses a list of dict nodes.
/// The input string is expected to contain the pattern name { ... }
/// \note The input string does not have the main enclosing parenthesis.
/// \note The string can contain spaces, commentaries and break lines.
/// \param s raw string containing the list
/// \return vector of points
template <>
inline std::vector<OpenFoamDict::DictNode>
OpenFoamDict::parseValuesFrom(const std::string &s) {
  HERMES_LOG("Parsing vector of dict nodes");
  std::vector<OpenFoamDict::DictNode> values;
  size_t i = 0;
  // read array size
  consumeAll(i, s);
  HERMES_ASSERT(std::isdigit(s[i]));
  auto count_token = nextToken(i, s);
  size_t count = std::stoull(count_token);
  if (count == 0)
    return values;
  auto nodes = readFields(i, s);
  for (const auto &node : nodes)
    values.emplace_back(node.second);
  if (values.size() < count)
    values.resize(count, values.empty() ? OpenFoamDict::DictNode{} : values[0]);
  HERMES_ASSERT(values.size() == count);
  return values;
}
/// \brief Parses a list of unsigned integers and stores in a vector of size_t.
/// The input string is expected to contain the pattern n v_1 ... v_n.
/// \note The input string does not have the main enclosing parenthesis.
/// \note The string can contain spaces, commentaries and break lines.
/// \param s raw string containing the list
/// \return vector of numbers
template <>
inline std::vector<size_t> OpenFoamDict::parseValuesFrom(const std::string &s) {
  HERMES_LOG("Parsing vector of size_t");
  size_t i = 0;
  std::vector<size_t> values;
  // read array size
  consumeAll(i, s);
  HERMES_ASSERT(std::isdigit(s[i]));
  auto count_token = nextToken(i, s);
  size_t count = std::stoull(count_token);
  if (count == 0)
    return values;
  // assume values are separated by lines
  while (i < s.size()) {
    consumeAll(i, s);
    auto v = nextToken(i, s);
    HERMES_ASSERT(hermes::Str::isInteger(v));
    values.emplace_back(std::stoull(v));
    consumeAll(i, s);
  }
  if (values.size() < count)
    values.resize(count, values.empty() ? 0 : values[0]);
  HERMES_ASSERT(values.size() == count);
  return values;
}
/// \brief Parses a list of integers and stores in a vector of integer.
/// The input string is expected to contain the pattern n v_1 ... v_n.
/// \note The input string does not have the main enclosing parenthesis.
/// \note The string can contain spaces, commentaries and break lines.
/// \param s raw string containing the list
/// \return vector of numbers
template <>
inline std::vector<i64> OpenFoamDict::parseValuesFrom(const std::string &s) {
  HERMES_LOG("Parsing vector of int");
  size_t i = 0;
  std::vector<i64> values;
  // read array size
  consumeAll(i, s);
  HERMES_ASSERT(std::isdigit(s[i]));
  auto count_token = nextToken(i, s);
  size_t count = std::stoull(count_token);
  if (count == 0)
    return values;
  // assume values are separated by lines
  while (i < s.size()) {
    consumeAll(i, s);
    auto v = nextToken(i, s);
    HERMES_ASSERT(hermes::Str::isInteger(v));
    values.emplace_back(std::stoi(v));
    consumeAll(i, s);
  }
  if (values.size() < count)
    values.resize(count, values.empty() ? 0 : values[0]);
  HERMES_ASSERT(values.size() == count);
  return values;
}
/// \brief Parses a list of unsigned integer lists and stores in a
/// vector<vector> of unsigned integers. The input string is expected to contain
/// the pattern n n_1(v v v) ... n_n(v v v).
/// \note The input string does not have the main enclosing parenthesis.
/// \note The string can contain spaces, commentaries and break lines.
/// \param s raw string containing the list
/// \return vector of lists
template <>
inline void
OpenFoamDict::parseValueListsFrom(const std::string &s,
                                  std::vector<std::vector<size_t>> &values) {
  HERMES_LOG("Parsing vector of faces");
  values.clear();
  size_t i = 0;
  // read array size
  consumeAll(i, s);
  HERMES_ASSERT(std::isdigit(s[i]));
  auto count_token = nextToken(i, s);
  size_t count = std::stoull(count_token);
  if (count == 0)
    return;
  while (i < s.size()) {
    // read a n(v v v v v) pattern
    consumeAll(i, s);
    // here we must find an integer
    auto t = nextToken(i, s);
    HERMES_ASSERT(hermes::Str::isInteger(t));
    auto n = std::stoul(t);
    consumeAll(i, s);
    // consume opening parenthesis
    HERMES_ASSERT(s[i] == '(');
    ++i;
    std::vector<size_t> vs;
    for (size_t j = 0; j < n; ++j) {
      consumeAll(i, s);
      t = nextToken(i, s);
      vs.emplace_back(std::stoull(t));
    }
    consumeAll(i, s);
    // consume closing parenthesis
    HERMES_ASSERT(s[i] == ')');
    ++i;
    values.emplace_back(vs);
    consumeAll(i, s);
  }
  if (values.size() < count)
    values.resize(count, values.empty() ? std::vector<size_t>{} : values[0]);
  HERMES_ASSERT(values.size() == count);
}

#endif // PSA_ANIM_VIS_OPENFOAM_VIEWER_OPENFOAM_PARSER_H
