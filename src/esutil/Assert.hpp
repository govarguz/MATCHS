#ifndef ASSERT_H
#define ASSERT_H

#include <string>

struct AssertError : public std::logic_error
{
    AssertError(std::string str)
      : std::logic_error(str)
    {

    }
};

// todo: after C++11, use std::to_string
#define ESPR_ASSERT(expression) if(!(expression)) throw AssertError(std::string("Assert in ") + __FILE__ + ":" + boost::lexical_cast<std::string>(__LINE__));


#endif // ASSERT_H
