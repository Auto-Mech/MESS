#include "key.hh"
#include <vector>
#include <stack>
#include <iomanip>

std::vector<std::vector<Key::_Val> > Key::_stack;

void Key::_check_stack ()
{
  Exception::Base funame = " Key::_check_stack (): ";

  if(!_stack.size())
    throw funame << "no key group has been initialized";
}

void Key::_check_range () const
{
  Exception::Base funame = " Key::_check_range (): ";

  _check_stack();

  if(_group >=_stack.size() || _index >= _stack[_group].size())
    throw funame << "out of range";
}

void Key::_add_key (const std::string& s)
{
  _check_stack();

  _group = _stack.size() - 1;
  _index = _stack.rbegin()->size();
  _stack.rbegin()->push_back(s);
}

bool Key::isinit () const
{
  _check_range();

  return _stack[_group][_index]._init;
}

bool Key::operator== (const std::string& s) const
{
  _check_range();

  if(s != _stack[_group][_index])
    return false;

  _stack[_group][_index]._init = true;
  return true;
}

bool Key::check_uninitialized_keys (std::ostream& to)
{
  _check_stack();

  bool res = true;
  for(std::vector<_Val>::const_iterator it = _stack.rbegin()->begin(); it != _stack.rbegin()->end(); ++it) {
    if(!it->_init) {
      if(res) {
	res = false;
	to << "The following parameters have not been initialized:";
      }
      to << "   " << *it;
    }
  }

  if(!res)
    to << "\n";

  return res;
}

void Key::show_all (std::ostream& to, int n)
{
  _check_stack();

  to << std::setw(n) << "" << "Available keys:\n";
  for(std::vector<_Val>::const_iterator it = _stack.rbegin()->begin();
      it != _stack.rbegin()->end(); ++it) {
    to << std::setw(n) << "" << *it << "\n";
  }
}

std::string Key::show_all (int n)
{
  std::ostringstream to;
  show_all(to, n);
  return to.str();
}
