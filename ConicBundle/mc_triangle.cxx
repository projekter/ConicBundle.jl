/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  mc_triangle.cxx
    This file is part of ConciBundle, a C/C++ library for convex optimization.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

***************************************************************************** */

#include <string.h>
#include <iostream>
#include <fstream>
#include <set>
#include <queue>
#include "CMsingleton.hxx"
#include "CMsymsparse.hxx"
#include "MatrixCBSolver.hxx"
#include "PSCAffineFunction.hxx"

using namespace CH_Matrix_Classes;
using namespace ConicBundle;
using namespace std;

//*****************************************************************************
//                                Triangle
//*****************************************************************************

/** @brief triangle inequalities are identified by indices i,j,k with
  0<=i<j<k<|nodes| and two signflafs signj,signk in {-1,1}. The
  inqualitiy reads
     signj*x_{ij}+signk*x_{ik}+signj*signk*x_{jk} >= -1. 
  For checking whether an inequality has been added already, this class
  stores the identifying triple and defines a sorting order. 
 */ 
class Triangle
{
private:
  Integer i,j,k; ///<  0<=i<abs(j)<abs(k) give the indices of the triangle inequality, the signs of j and k give the type of the triangle inequality
  
public:
  
  ///construct
  Triangle(Integer ii,Integer jj,Integer kk,Integer signj,Integer signk):i(ii),j(jj*signj),k(kk*signk)
  { assert((0<=ii)&&(ii<jj)&&(jj<kk)&&(std::abs(signj)==1)&&(std::abs(signk)==1)); }

  ///copy constructor 
  Triangle(const Triangle& tr)
  { i=tr.i; j=tr.j; k=tr.k; }

  ///destructor (nothing to do)
  ~Triangle(){}

  ///comparison operator for standard template library
  bool operator<(const Triangle& tr) const
  {
    if (i!=tr.i) return (i<tr.i);
    if (std::abs(j)!=std::abs(tr.j)) return (std::abs(j)<std::abs(tr.j));
    if (std::abs(k)!=std::abs(tr.k)) return (std::abs(k)<std::abs(tr.k));
    if (j!=tr.j) return (j<tr.j);
    return (k<tr.k);
  }

  ///retrieve triangle data
  void get_ijk(Integer& ii,Integer& jj,Integer& kk,Integer& signj,Integer& signk) const
  { ii=i; jj=std::abs(j); kk=std::abs(k); signj=(j<0?-1:1); signk=(k<0?-1:1); }
};

//*****************************************************************************
//                                TriangleViolation
//*****************************************************************************

/** @brief in order to find the most violated triangles by a priority queue this stores a triangle together with its violation and allows sorting according to the violation (-1. is still ok, the more negative the more violated)
 */ 
class TriangleViolation
{
private:
  /// the triangle
  Triangle triangle;
  /// the violation of the triangle inequality, more negative is more violated
  Real violation;
public:
  TriangleViolation(Integer ii,Integer jj,Integer kk,Integer signj,Integer signk,Real viol):
    triangle(ii,jj,kk,signj,signk),violation(viol)
  {}

  TriangleViolation(const TriangleViolation& Tr):
    triangle(Tr.triangle),violation(Tr.violation)
  {}

  bool operator<(const TriangleViolation& Tr) const
  {
    return (violation<Tr.violation);
  }

  const Triangle& get_triangle() const
  { return triangle; }

  Real get_violation() const
  {return violation;}
};


//*****************************************************************************
//                               added_triangles
//*****************************************************************************

//for storing information about the added triangle inequalities 
typedef std::set<Triangle> AddedTriangles;


//*****************************************************************************
//                      update_triangle_constraints()
//*****************************************************************************

/** @brief retrieves a primal approximation X and finds the most violated 
  triangle ineuqalities for this. If they are not yet include, i.e., if 
  they are not yet in added_triangles, they are added to the problem,
  i.e., to the solver and via the solver to the function. 

  Note, it is much more efficient to separate odd cycles with respect 
  to a dynamically extended sparse support, but for an example this would 
  be too complicated to show.
*/

int update_triangle_constraints(MatrixCBSolver& solver,PSCAffineFunction& mc,AddedTriangles& added_triangles,const Indexmatrix* fixed_indices,Integer nnodes,const GramSparsePSCPrimal* primalX,Real& violation)
{
  //-- get and form an approximation of the primal matrix X
  Symmatrix X(*primalX);
  rankadd(primalX->get_grammatrix(),X,1.,1.);

  //-- find the n most violated triangle inequalities
  std::size_t at_most_n_violated=std::size_t(max(nnodes,10));
  Real violation_threshold=-1-1e-3;
  
  std::priority_queue<TriangleViolation> most_violated;
  for(Integer i=0;i<nnodes;i++){
    for(Integer j=i+1;j<nnodes;j++){
      Real xij=X(i,j);
      for(Integer k=j+1;k<nnodes;k++){
	Real xik=X(i,k);
	Real xjk=X(j,k);
	int signj=1;
	int signk=1;
	Real viol=xij+xik+xjk;
	Real v=xij-xik-xjk;
	if (v<viol){
	  viol=v;
	  signk=-1;
	}
	v=-xij+xik-xjk;
	if (v<viol){
	  viol=v;
	  signj=-1;
	  signk=1;
	}
	v=-xij-xik+xjk;
	if (v<viol){
	  viol=v;
	  signj=-1;
	  signk=-1;
	}
	if (viol<violation_threshold){
	  if (most_violated.size()==at_most_n_violated){
	    most_violated.pop();
	  }
	  most_violated.push(TriangleViolation(i,j,k,signj,signk,viol));
	  if (most_violated.size()==at_most_n_violated){
	    violation_threshold=most_violated.top().get_violation();
	  }
	}
      }  // end k
    } // end j
  } // end i

  //-- build the coefficient matrices for those that are not yet included
  Indexmatrix indi(3,1,Integer(0));  //row indices for sparse representation
  Indexmatrix indj(3,1,Integer(0));  //col indices for sparse representation
  Matrix val(3,1,1.);                //values for sparse representation
  Sparsesym A;                       //sparse symmetric matrix representation
  CoeffmatVector tr_ineqs;           //collects the ineqs coefficient matrices
  tr_ineqs.reserve(most_violated.size());
  Real sumviol=0.;     //for computing the average violation of added ineqs
  
  while(most_violated.size()>0){
    //check whether the inequality is new
    if (added_triangles.find(most_violated.top().get_triangle())==added_triangles.end()){
      //it is new, form and append it multiplied by -1 for nonneg. multipiers
      violation=-1.-most_violated.top().get_violation(); //the last is the most violated one
      sumviol+=violation;
      added_triangles.insert(most_violated.top().get_triangle());
      Integer i,j,k,signj,signk;
      most_violated.top().get_triangle().get_ijk(i,j,k,signj,signk);
      indi(0)=i; indj(0)=j; val(0)=.5*signj;
      indi(1)=i; indj(1)=k; val(1)=.5*signk;
      indi(2)=j; indj(2)=k; val(2)=.5*signj*signk;
      Sparsesym A(nnodes,3,indi,indj,val);
      tr_ineqs.push_back(new CMsymsparse(A));
    }
    most_violated.pop();
  }

  //-- add the coefficient matrices (the primal constraints) to the problem
  Integer ncols(Integer(tr_ineqs.size()));
  if (ncols>0){
    std::cout<<" nnew="<<ncols<<" maxviol="<<violation<<" avgviol="<<sumviol/ncols<<std::endl;
    //form the row of the coefficient matrix that has to be appended
    indi.init(ncols,1,Integer(0));
    indj.init(Range(0,ncols-1));
    SparseCoeffmatMatrix scm(mc.get_opAt().blockdim(),ncols,&indi,&indj,&tr_ineqs);
    //form the modification info for the PSCAffineFunction mc
    PSCAffineModification funmod(solver.get_dim(),mc.get_opAt().blockdim(),&solver); 
    funmod.add_append_vars(ncols,&scm);
    funmod.set_skip_extension(true); //constraints likely not in support of the Laplacian
    //set the function-object-modifcation-map
    FunObjModMap modmap;
    modmap[&mc]=&funmod;
    //append the new Lagrange mulitpliers for <=1. constraints 
    Matrix lower_bounds(ncols,1,0.);
    Matrix costs(ncols,1,1.);
    if (solver.append_variables(ncols,&lower_bounds,0,0,0,&costs,&modmap)){
      cout<<"**** ERROR in solver.append_variables(...)"<<endl;
      return 1;
    }
  }

  //---- delete the triangle constraints whose multipliers are fixed to zero
  Matrix centery;
  if (solver.get_center(centery)){
    cout<<"**** ERROR in solver.get_candidate()"<<endl;
    return 1;
  }
  if (fixed_indices) {
    Indexmatrix delind(fixed_indices->rowdim(),1);
    delind.init(0,1,Integer(0));
    for(Integer i=nnodes;i<fixed_indices->rowdim();i++){
      if ((*fixed_indices)(i)&&(centery(i)==0.)){
	delind.concat_below(i);
      }
    }

    if (delind.rowdim()>0){
      //check which are triangle constraints for cleaning up AddedTriangles
      Indexmatrix indi,indj;
      Matrix val;
      for(Integer i=0;i<delind.rowdim();i++){
	const CMsymsparse* cm=dynamic_cast<const CMsymsparse*>(mc.get_opAt()(0,delind(i)).ptr());
	assert(cm);
	cm->sparse(indi,indj,val,2.);
	assert(indi.rowdim()==3);
	Integer ii=indi(0);
	Integer jj=indj(0);
	Integer kk=indj(1);
	assert((indi(1)==ii)&&(indi(2)==jj)&&(indj(2)==kk)&&(ii<jj)&&(jj<kk));
	Integer signj=val(0)>0?1:-1;
	Integer signk=val(1)>0?1:-1;
	assert(std::fabs(val(0)-signj)<1e-6);
	assert(std::fabs(val(1)-signk)<1e-6);
	assert(std::fabs(val(2)-signj*signk)<1e-6);
	AddedTriangles::iterator it=added_triangles.find(Triangle(ii,jj,kk,signj,signk));
	assert(it!=added_triangles.end());
	added_triangles.erase(it);
      }
      
      //delete the variable
      Indexmatrix map_to_old;
      if(solver.delete_variables(delind,map_to_old)){
	cout<<"**** ERROR in solver.get_lbounds()"<<endl;
	return 1;
      }
      
      cout<<" deleted "<<delind.rowdim();
    }
  }
  cout<<" next_dim="<<solver.get_dim()<<std::endl;

  return 0;
}


//*****************************************************************************
//                               GW_rounding()
//*****************************************************************************

/** @brief use Goemans-Williamson random hyperplane rounding on the approximate primal solution to generate a random cut and store the best one found.
 
   No local improvement heuristic is employed here to keep the code simple.
*/

int GW_rounding(const GramSparsePSCPrimal* primalX,const Sparsesym& L,Matrix& bestcut,Real& bestcut_val)
{
  const Matrix& Gram_mat=primalX->get_grammatrix();
  //--- try the signs of the most important column first
  Matrix cut(Gram_mat.col(0));
  for(Integer i=0;i<cut.rowdim();i++)
    cut(i)=cut(i)>0.?1:-1;
  Real cutval=ip(cut,L*cut);
  if (cutval>bestcut_val){
    bestcut=cut;
    bestcut_val=cutval;
  }
  //--- try Goemans Williamson rounding
  Matrix rand_dir;
  Integer n_tries=max(L.rowdim()/10,10);
  for (Integer r=0;r<n_tries;r++){
    rand_dir.rand_normal(Gram_mat.coldim(),1);
    genmult(Gram_mat,rand_dir,cut);
    for(Integer i=0;i<cut.rowdim();i++)
      cut(i)=cut(i)>0.?1:-1;
    cutval=ip(cut,L*cut);
    if (cutval>bestcut_val){
      bestcut=cut;
      bestcut_val=cutval;
    }
  }
  cout<<" best_cut="<<bestcut_val<<std::endl;
  return 0;
}

//*****************************************************************************
//                                 main()
//*****************************************************************************

int main(int argc,char **argv)
{
  //---- open the input file or get the input from stdin
  ifstream fin;
  istream* in=&cin;
  if (argc!=2){
    cerr<<" usage: mc_triangle <filename>\n";
    cerr<<" for reading from stdin use: mc_triangle - \n";
    cerr<<" input format [node indices start at 1]:  nnodes medges edge1_node1 edge1_node2 edge1_val ... edgem_node1 edgem_node2 edgem_val"<<endl;
    return 1;
  }
  if (strcmp("-",argv[1])!=0){
    //open filename
    fin.open(argv[1]);
    if (!fin.good()){
      cerr<<" failure in opening file named "<<argv[1]<<endl;
      return 1;
    }
    in=&fin;
  }
  
  //---- read the graph into (indi,indj,val) triples (replace val by 1.)
  Integer nnodes,medges;
  
  (*in)>>nnodes>>medges;
  
  Indexmatrix indi(medges,1,Integer(0));
  Indexmatrix indj(medges,1,Integer(0));
  Matrix val(medges,1,1.);
  
  for (int i = 0; i <medges; i++) {
    int head;
    int tail;
    double d;
    (*in)>>tail>>head>>d;      
    indi(i)=tail-1;
    indj(i)=head-1;
    //val(i)=d;               //edge weight d could be plugged in here
  }
  fin.close();
  
  //---- form Laplacian/4. as a Sparse Symmetric cost matrix
  Sparsesym L(nnodes,medges,indi,indj,val);
  Matrix Ldiag(diag(L));
  L-=sparseDiag(Ldiag);
  Ldiag.init(nnodes,1,sum(L)/Real(nnodes));
  L*=-1;
  L+=sparseDiag(Ldiag);
  L/=4.;
  Matrix bestcut(nnodes,1,1.);                //initialize with empty cut
  Real bestcut_val=ip(bestcut,L*bestcut);

  //---- form the Affine Matrix Function Oracle
  Indexmatrix Xdim(1,1,nnodes);
  SparseCoeffmatMatrix C(Xdim,1);
  C.set(0,0,new CMsymsparse(L));

  SparseCoeffmatMatrix opAt(Xdim,nnodes);
  for (Integer i=0;i<nnodes;i++)
    opAt.set(0,i,new CMsingleton(nnodes,i,i,-1.));
  
  PSCAffineFunction mc(C,opAt,new GramSparsePSCPrimal(L));
  mc.set_out(&cout,0);

  //---- initialize the solver and the problem
  MatrixCBSolver cbsolver(&cout,1);
  Matrix lb(nnodes,1,CB_minus_infinity);
  Matrix ub(nnodes,1,CB_plus_infinity);
  Matrix rhs(nnodes,1,1.);
  cbsolver.init_problem(nnodes,&lb,&ub,0,&rhs);
  cbsolver.add_function(mc,Real(nnodes),ObjectiveFunction,0,true); 
   
  //---- call the solver
  AddedTriangles added_triangles;  //keeps track of added triangle inequalities
  Real violation=0.;               //latest maximum violation of a tr.-ineq.
  // if dual variables to tr.-ineqs. stay zero allow to fix them there 
  cbsolver.set_active_bounds_fixing(true);

  do {
    
    // solve with option to stop after 10 null steps or the next descent step
    int retval=cbsolver.solve(10,true);
    if (retval) {
      cout<<"**** ERROR cbsolver.solve(..) returned "<<retval; 
    }

    //retrieve the current approximat primal semidefinite solution
    const GramSparsePSCPrimal* primalX= dynamic_cast<const GramSparsePSCPrimal*>(cbsolver.get_approximate_primal(mc));
    if (primalX==0){
      cout<<"**** ERROR in cbsolver.get_approximate_primal()"<<endl;
      return 1;
    }

    //try to improve the best found solution by Goemans-Williamson rounding
    retval=GW_rounding(primalX,L,bestcut,bestcut_val);
    if (retval) {
      cout<<"**** GW_rounding(..) returned "<<retval; 
    }
    if (cbsolver.get_objval()-bestcut_val<1.-1e-6){
      cout<<" gap between bound and feasible solution <1."<<std::endl;
      break;
    }
    
    //retrieve incidence vector of fixed variables
    const Indexmatrix* fixed_indices=cbsolver.get_fixed_active_bounds();

    // add new triangle inequalities and delete fixed ones
    violation=0;
    retval=update_triangle_constraints(cbsolver,mc,added_triangles,fixed_indices,nnodes,primalX,violation);
    if (retval) {
      cout<<"**** ERROR spearate(....) returned "<<retval; 
    }
    
  } while((cbsolver.termination_code()==0)||
	  ((cbsolver.termination_code()<2)&&(violation>1e-2)));

  //---- print some information about the solution (process) and the problem
  cbsolver.print_termination_code(cout);

  return 0;
}


