/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#include "logical.hh"
#include "io.hh"

#include <list>
#include <cctype>

bool Logical::BinExpr::evaluate (const std::vector<bool>& l) const  throw(Error::General)
{
  const char funame [] = "Logical::BinExpr::evaluate: ";
    
  switch(_op) {
  case OR: 
    return _x1->evaluate(l) || _x2->evaluate(l);
  case AND:
    return _x1->evaluate(l) && _x2->evaluate(l);
  default:
    std::cerr << funame << "unknown operation " << _op << "\n";
    throw Error::Range();
  }
}

void Logical::VarExpr::init (const std::map<std::string, int>& l) throw(Error::General)
{
  const char funame [] = "Logical::VarExpr::init: ";

  typedef std::map<std::string, int>::const_iterator Mit;

  Mit mit = l.find(_name);

  if(mit == l.end()) {
    std::cerr << funame << "no name <" << _name << "> in the list: ";
    for(mit = l.begin(); mit != l.end(); ++mit)
      std::cerr << mit->first << " ";
    std::cerr << "\n";
    throw Error::Init();
  }

  _var = mit->second;
  if(_var < 0) {
    std::cerr << funame << "index should not be negative\n";
    throw Error::Init();
  }
}

SharedPointer<Logical::Expr> Logical::read_expr (std::istream& from) throw(Error::General)
{
    char next;
    
    std::list<SharedPointer<Expr> > x;
    std::list<BinExpr::Oper> op;
    while(1) {
	x.push_back(read_term(from));
	next = IO::skip_space(from);
	if(next == '|') 
	  op.push_back(BinExpr::OR);
	else if(next == '&')
	  op.push_back(BinExpr::AND);
	else
	  break;
	from.get(next);
    }

    // collapse all AND operations
    std::list<SharedPointer<Expr> >::iterator  xit, xit1;
    xit= x.begin(); ++xit;
    std::list<BinExpr::Oper>::iterator opit = op.begin();

    while(opit != op.end())
      if(*opit == BinExpr::AND) {
	xit1 = xit;
	--xit1;
	*xit1 = *xit1 & *xit;
	xit  =  x.erase(xit);
	opit = op.erase(opit);
      }
      else {
	++xit;
	++opit;
      }	
    
    xit = x.begin();
    SharedPointer<Expr> res = *xit;
    while(++xit != x.end())
      res = res | *xit;
    return res;
}


SharedPointer<Logical::Expr> Logical::read_term (std::istream& from) throw(Error::General)
{
  const char funame [] = "Logical::read_term: ";

  SharedPointer<Expr> res;
  char next = IO::skip_space(from);

  if(!from) {
    std::cerr << funame << "stream is corrupted\n";
    throw Error::Form();
  }

  switch(next) {
  case '!':
    from.get(next); 
    return negate(read_term(from));
  case '(':
    from.get(next);
    res = read_expr(from);
    from.get(next);
    if(next != ')') {
      std::cerr << funame << "no closing parenthesis\n";
      throw Error::Form();
    }
    return res;
  default:
    return read_name(from);
  }
}

SharedPointer<Logical::Expr> Logical::read_name (std::istream& from) throw(Error::General)
{
  const char funame [] = "Logical::read_name: ";

  std::string name; 
  char next = from.get();

  if(!std::isalpha(next)) {
    std::cerr << funame << "should begin with letter\n";
    throw Error::Form();
  }

  name = next;
  while(from.get(next) && std::isalnum(next))
    name += next;
  
  if(from)
    from.putback(next);
      
  return SharedPointer<Expr>(new VarExpr(name));
}
