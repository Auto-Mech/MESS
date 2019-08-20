/*
        Chemical Kinetics and Dynamics Library
        Copyright (C) 2008-2013, Yuri Georgievski <ygeorgi@anl.gov>

        This library is free software; you can redistribute it and/or
        modify it under the terms of the GNU Library General Public
        License as published by the Free Software Foundation; either
        version 2 of the License, or (at your option) any later version.

        This library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
        Library General Public License for more details.
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
       = 0;
    virtual void init (const std::map<std::string, int>&)
       = 0;
  };

  SharedPointer<Expr> read_expr (std::istream& from) ;
  SharedPointer<Expr> read_term (std::istream& from) ;
  SharedPointer<Expr> read_name (std::istream& from) ;

  // unary negation expression 
  class UniExpr : public Expr
  {
    SharedPointer<Expr> _x;

    UniExpr (SharedPointer<Expr> x) : _x(x) {}
  public:
 
    bool evaluate (const std::vector<bool>& l) const 
       { return !(_x->evaluate(l)); }

    void init (const std::map<std::string, int>& l)
       { _x->init(l); }

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
      ;

    void init (const std::map<std::string, int>& l)
       { _x1->init(l); _x2->init(l); }

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
      ;
    void init (const std::map<std::string, int>&)
      ;
  };

  inline bool VarExpr::evaluate (const std::vector<bool>& l) const 
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
