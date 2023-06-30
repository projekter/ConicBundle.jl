/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/ModificationTreeData.hxx
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




#ifndef CONICBUNDLE_MODIFICATIONTREEDATA_HXX
#define CONICBUNDLE_MODIFICATIONTREEDATA_HXX


/**  @file ModificationTreeData.hxx
    @brief Header declaring the class ConicBundle::ModificationTreeData
    @version 1.0
    @date 2014-11-10
    @author Christoph Helmberg
*/

#include "AffineFunctionTransformation.hxx"
#include "FunctionObjectModification.hxx"

namespace ConicBundle {

/** @ingroup dynamic_modification_support
*/
//@{

class SumBlockModel;
class ModificationTreeData;

/** @brief this serves either as searchable storage for nodes on the heap or as an adjacency list of adjacent nodes  
 */  
typedef std::map<const FunctionObject*, ModificationTreeData*> FunctionMap;


/** @brief Represents a tree node of the tree within MatrixCBSolver that describes the problem with its interacting functions and that allows to collect modifications before executing them on the respective models

 */

class ModificationTreeData : public CBout
{    
private:
  //-----------------------------------------------------------------------
  /** @name  data for describing function properties and its place in the tree
   */
  //@{

  const FunctionObject* funobject; ///< the actual function object, i.e. the oracle
  FunctionObject* wrapper; ///< if !=funobject, the function is wrapped in *wrapper as an oracle interface 
  SumBlockModel* model; ///< the model used for funobject or its wrapper
  
  ModificationTreeData* parent;   ///< if the function is part of parent function (sum,...), this points to it
  FunctionMap children; ///< if the function consists of several subfunctions, these are pointed to here

  //@}
    
  //-----------------------------------------------------------------------
  /** @name  data needed for modifications
   */
  //@{
  
  CH_Matrix_Classes::Integer fixed_dimension; ///< if negative, this means not fixed, otherwise this is the number of parameters and they may never be increased, reordered or reassigned 
  AFTModification aftmod; ///< this also gives an account of all modifications of the ground set of the oracle via its rows 
  OracleModification* oraclemod;
  
  //@}

  //-----------------------------------------------------------------------

public:

  /** @brief initialize the function object with its dimension, possibly with AFT and dynamic settings, its model, possibly with an oracle wrapper if required, and output settings
      
      
   */
   ModificationTreeData(const FunctionObject* fo,
			  FunctionObject* wr,
			  SumBlockModel* sbm,
			  CH_Matrix_Classes::Integer groundset_dim,
			  CH_Matrix_Classes::Integer fun_dim_fixed,
			  AffineFunctionTransformation* aft,
			  const CBout* cb);

  /// destructor, does not delete the oracle but the model and the wrapper (if there is one) and assumes that all children have been removed already
  ~ModificationTreeData();

  /// returns the model
  SumBlockModel* get_model() const {return model;}

  /// returns true if some modifications have been collected and need to be executed
  bool modification_pending() const;
    
  /// returns true if some modifications have been collected and need to be executed
  bool 
  pending_oracle_modification(
			      CH_Matrix_Classes::Integer& old_dim,
			      CH_Matrix_Classes::Integer& new_dim,
			      CH_Matrix_Classes::Integer& append_dim,
			      const CH_Matrix_Classes::Indexmatrix*& map_to_old,
			      const CH_Matrix_Classes::Indexmatrix*& deleted_indices,
			      const CH_Matrix_Classes::Indexmatrix*& new_indices,
			      const OracleModification*& oracle_modification
			      ) const;
    
  /// add a child (possibly holding a subtree) to this node and to the corresponding model; if any error occurs the node/subtree is not added
  int add_child(ModificationTreeData* fmd);
    
  /// remove the subtree rooted at this node from the parent of this but do not delete it
  int unlink_subtree();
        
  /// recursively remove all tree nodes of the subtree (except for this) from funmap and destruct/delete them 
  int delete_descendants(FunctionMap& funmap);
    
  /// Recursively add the parent node of this to funmodmap if it is not yet in there (the purpose is to propagate necessary modifications to the root)
  int add_parents_to_map(FunObjModMap& funmodmap);
    
  /** @brief starting from the root recursively collect the required modifications in each node and propagate them to the children; 

      if modmap specifies modifications for this function object, they
      override add_dim and map_to_old. If no special modifications are given,
      add_dim is carried out before map_to_old. If add_dim==0 or map_to_old==0
      there are no modifications of htis kind.

      Independent of whether modifications at this level exist, 
      the children are always called afterwards with the modifications
      resulting from the step at this level.
  */ 
  int update_subtree_modification(CH_Matrix_Classes::Integer add_dim,
				  const CH_Matrix_Classes::Indexmatrix* map_to_old,
				  const FunObjModMap* modmap);

  /// recursively add modifications whenever they deviate from default modifications referenced by wrapper instead of funobject whenever it exists 
  int collect_subtree_modification(FunObjModMap& modmap);

  /// set the data to the state achieved by executing the stored modifications, no modifications are pending afterwards (this does not change the models but only the modification data)
  int subtree_modification_performed();

  /// recursively clear all pending modifications without being executed, i.e., go back to the state without modifications 
int clear_subtree_modification();

  /// set output levels for this, wrapper, model and aftmod
  void set_cbout(const CBout* cb,int incr=-1);

};


  //@}

}

#endif

