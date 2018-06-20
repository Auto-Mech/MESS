/*
    Copyright (C) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
    Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

    See https://github.com/PACChem/MESS for copyright and licensing details.
*/

#ifndef LOGICAL_HH
#define LOGICAL_HH

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include "shared.hh"

namespace Logical 
{
  class Expr 
  {
  public:
    virtual bool evaluate (const std::vector<bool>&) const 
      throw(Error::General) = 0;
    virtual void init (const std::map<std::string, int>&)
      throw(Error::General) = 0;
  };

  SharedPointer<Expr> read_expr (std::istream& from) throw(Error::General);
  SharedPointer<Expr> read_term (std::istream& from) throw(Error::General);
  SharedPointer<Expr> read_name (std::istream& from) throw(Error::General);

  // unary negation expression 
  class UniExpr : public Expr
  {
    SharedPointer<Expr> _x;

    UniExpr (SharedPointer<Expr> x) : _x(x) {}
  public:
 
    bool evaluate (const std::vector<bool>& l) const 
      throw(Error::General) { return !(_x->evaluate(l)); }

    void init (const std::map<std::string, int>& l)
      throw(Error::General) { _x->init(l); }

    friend SharedPointer<Expr> negate (SharedPointer<Expr>);
  };

  inline SharedPointer<Expr> negate (SharedPointer<Expr> x) 
  { return SharedPointer<Expr>(new UniExpr(x)); }

  // binary expression
  class BinExpr : public Expr 
  {
  public:
    enum Oper {OR, AND};

  private:
    Oper _op;
    SharedPointer<Expr> _x1;
    SharedPointer<Expr> _x2;

    BinExpr (SharedPointer<Expr> x1, SharedPointer<Expr> x2, Oper op) : _x1(x1), _x2(x2) , _op(op) {}

  public:
    bool evaluate (const std::vector<bool>&) const 
      throw(Error::General);

    void init (const std::map<std::string, int>& l)
      throw(Error::General) { _x1->init(l); _x2->init(l); }

    friend SharedPointer<Expr> operator& (SharedPointer<Expr>, SharedPointer<Expr>); 
    friend SharedPointer<Expr> operator| (SharedPointer<Expr>, SharedPointer<Expr>); 
  };

  inline SharedPointer<Expr> operator& (SharedPointer<Expr> x1, SharedPointer<Expr> x2) 
  { return SharedPointer<Expr>(new BinExpr(x1, x2, BinExpr::AND)); }

  inline SharedPointer<Expr> operator| (SharedPointer<Expr> x1, SharedPointer<Expr> x2) 
  { return SharedPointer<Expr>(new BinExpr(x1, x2, BinExpr::OR)); }

  // variable expression
  class VarExpr : public Expr
  {
    std::string _name;
    int _var;
    
  public:
    VarExpr(const std::string& s) : _name(s), _var(-1) {}

    bool evaluate (const std::vector<bool>&) const
      throw(Error::General);
    void init (const std::map<std::string, int>&)
      throw(Error::General);
  };

  inline bool VarExpr::evaluate (const std::vector<bool>& l) const throw(Error::General)
  {
    const char funame [] = "Logical::VarExpr::evaluate: ";

#ifdef DEBUG
    if(_var < 0) {
      std::cerr << funame << "not initialized\n";
      throw Error::Init();
    }
    if(_var >= l.size()) {
      std::cerr << funame << "out of range\n";
      throw Error::Range();
    }
#endif

    return l[_var]; 
  }

}// Logical

#endif
