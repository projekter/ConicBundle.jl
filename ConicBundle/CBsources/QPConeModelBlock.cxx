/* ****************************************************************************

    Copyright (C) 2004-2021  Christoph Helmberg

    ConicBundle, Version 1.a.2
    File:  CBsources/QPConeModelBlock.cxx
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


#include "QPConeModelBlock.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

  Real QPConeModelBlock::evaluate_trace(Matrix& vec) {
    return ip(vec, trace_vec);
  }

  void QPConeModelBlock::modelx_changed() {
    QPModelBlock::modelx_changed();
    sysrhs_model.init(0, 1, 0.);
    sysrhs_trace = 0.;
    sysinv_trace.init(0, 1, 0.);
    Btsysinv_trace.init(0, 1, 0.);
    sys_trace = -1.;
    schur_trace = -1.;
    complrhs_trace = 0;
  }

  void QPConeModelBlock::clear() {
    modeldim = 0;
    mu_dim = 0;
    ft = ObjectiveFunction;
    trace_vec.init(0, 1, 0.);
    last_rhs_mu = 0.;
    mu = 0.;
    old_mu = 0.;
    trace_rhs = 0.;
    trace_slack = 0;
    trace_dual = 0;
    trace_delta_slack = 0.;
    trace_delta_dual = 0.;
    sys_trace = -1.;
    complrhs_trace = 0.;
    diff_trace = 0.;
    sysrhs_model.init(0, 1, 0.);
    sysrhs_trace = 0;
    sysinv_trace.init(0, 1, 0.);
    for (unsigned i = 0; i < block.size(); i++) {
      delete block[i];
    }
    block.clear();
    nncblock = 0;
    boxblock = 0;
    socblock.clear();
    pscblock.clear();
    QPModelBlock::clear();

    oracle_data = 0;
  }

  QPConeModelBlock::QPConeModelBlock(CBout* cb, int cbinc) :QPModelBlock(cb, cbinc), QPConeModelDataObject(cb, cbinc), nncblock(0), boxblock(0) {
    clear();
  }

  QPConeModelBlock::~QPConeModelBlock() {
    clear();
  }

  void QPConeModelBlock::recursive_delete_and_clear() {
    clear();
  }

  int QPConeModelBlock::recursive_copy_data_of(QPModelBlockObject* p) {
    QPConeModelBlock* cp = dynamic_cast<QPConeModelBlock*>(p);
    if ((cp == 0) || (cp->block.size() != block.size()))
      return 1;

    constant_minorant = cp->constant_minorant;
    bundle = cp->bundle;
    modelx = cp->modelx;
    Bt = cp->Bt;
    modeldx = cp->modeldx;
    modeldcstr = cp->modeldcstr;
    sysviol_model = cp->sysviol_model;
    sysviol_constraints = cp->sysviol_constraints;
    modelx_aggregate = cp->modelx_aggregate;

    modeldim = cp->modeldim;
    mu_dim = cp->mu_dim;
    ft = cp->ft;
    use_trace = cp->use_trace;
    last_rhs_mu = cp->last_rhs_mu;
    mu = cp->mu;
    old_mu = cp->old_mu;
    trace_vec = cp->trace_vec;
    trace_rhs = cp->trace_rhs;
    trace_slack = cp->trace_slack;
    trace_dual = cp->trace_dual;
    trace_delta_slack = cp->trace_delta_slack;
    trace_delta_dual = cp->trace_delta_dual;
    sys_trace = cp->sys_trace;
    complrhs_trace = cp->complrhs_trace;
    diff_model = cp->diff_model;
    diff_trace = cp->diff_trace;
    sysrhs_model = cp->sysrhs_model;
    sysrhs_trace = cp->sysrhs_trace;
    sysinv_trace = cp->sysinv_trace;
    Btsysinv_trace = cp->Btsysinv_trace;
    schur_trace = cp->schur_trace;

    int err = 0;
    for (unsigned i = 0; i < block.size(); i++) {
      err += block[i]->copy_from(cp->block[i]);
    }

    oracle_data = cp->oracle_data;

    return err;
  }

  /// return a cloned object on the heap
  QPModelBlockObject* QPConeModelBlock::clone() {
    QPConeModelBlock* p = new QPConeModelBlock(this, 0);

    p->constant_minorant = constant_minorant;
    p->bundle = bundle;
    p->modelx = modelx;
    p->Bt = Bt;
    p->modeldx = modeldx;
    p->modeldcstr = modeldcstr;
    p->sysviol_model = sysviol_model;
    p->sysviol_constraints = sysviol_constraints;
    p->modelx_aggregate = modelx_aggregate;

    p->modeldim = modeldim;
    p->mu_dim = mu_dim;
    p->ft = ft;
    p->use_trace = use_trace;
    p->last_rhs_mu = last_rhs_mu;
    p->mu = mu;
    p->old_mu = old_mu;
    p->trace_vec = trace_vec;
    p->trace_rhs = trace_rhs;
    p->trace_slack = trace_slack;
    p->trace_dual = trace_dual;
    p->trace_delta_slack = trace_delta_slack;
    p->trace_delta_dual = trace_delta_dual;
    p->sys_trace = sys_trace;
    p->complrhs_trace = complrhs_trace;
    p->diff_model = diff_model;
    p->diff_trace = diff_trace;
    p->sysrhs_model = sysrhs_model;
    p->sysrhs_trace = sysrhs_trace;
    p->sysinv_trace = sysinv_trace;
    p->Btsysinv_trace = Btsysinv_trace;
    p->schur_trace = schur_trace;

    p->nncblock = 0;
    if (nncblock) {
      p->nncblock = dynamic_cast<NNCIPBundleBlock*>(nncblock->clone());
      p->block.push_back(p->nncblock);
    }
    for (unsigned i = 0; i < socblock.size(); i++) {
      SOCIPBundleBlock* ip = dynamic_cast<SOCIPBundleBlock*>(socblock[i]->clone());
      p->socblock.push_back(ip);
      p->block.push_back(ip);
    }
    for (unsigned i = 0; i < pscblock.size(); i++) {
      PSCIPBundleBlock* ip = dynamic_cast<PSCIPBundleBlock*>(pscblock[i]->clone());
      p->pscblock.push_back(ip);
      p->block.push_back(ip);
    }
    p->boxblock = 0;
    if (boxblock) {
      p->boxblock = dynamic_cast<BoxIPBundleBlock*>(boxblock->clone());
      p->block.push_back(p->boxblock);
    }

    p->oracle_data = oracle_data;

    return p;
  }



  //------------ QPConeModelBlockObject routines

  int QPConeModelBlock::init(const MinorantPointer& in_constant_minorant,
    const MinorantBundle& in_bundle,
    Integer nnc_dim,
    const Indexmatrix& soc_dim,
    const Indexmatrix& psc_dim,
    const Matrix& box_lb,
    const Matrix& box_ub,
    Real b,
    FunctionTask inft,
    QPModelOracleDataObject* in_oracle_data,
    bool scale_box) {
    clear();

    constant_minorant.push_back(in_constant_minorant);
    bundle.push_back(in_bundle);


    trace_rhs = b;
    ft = inft;
    oracle_data = in_oracle_data;

    // TEST output begin
    // Integer bsz=Integer(bundle.size());
    // Matrix tmpvec(bsz,1,0.);
    // Matrix tmpmat(2,bsz,0.);
    // for(Integer i=0;i<bsz;i++){
    //   bundle[unsigned(i)].get_minorant(tmpvec(i),tmpmat,i);
    // }
    // std::cout<<" B0="<<tmpvec<<std::endl;
    // std::cout<<" B="<<tmpmat<<std::endl;    
    // bsz=Integer(in_bundle.size());
    // tmpvec.init(bsz,1,0.);
    // tmpmat.init(2,bsz,0.);
    // for(Integer i=0;i<bsz;i++){
    //   in_bundle[unsigned(i)].get_minorant(tmpvec(i),tmpmat,i);
    // }
    // std::cout<<" inB0="<<tmpvec<<std::endl;
    // std::cout<<" inB="<<tmpmat<<std::endl;    
    // TEST output end

    modeldim = 0;
    mu_dim = 0;
    assert(nnc_dim >= 0);
    use_trace = false;
    if (nnc_dim > 0) {
      use_trace = true;
      modeldim += nnc_dim;
      mu_dim += nnc_dim;
      nncblock = new NNCIPBundleBlock(nnc_dim, this);
      block.push_back(nncblock);
      nncblock->set_cbout(this);
      nncblock->set_oracle_data(oracle_data);
    }
    for (Integer i = 0; i < soc_dim.dim(); i++) {
      assert(soc_dim(i) > 1);
      use_trace = true;
      modeldim += soc_dim(i);
      mu_dim += 1;
      SOCIPBundleBlock* sb = new SOCIPBundleBlock(soc_dim(i), this);
      socblock.push_back(sb);
      block.push_back(sb);
      sb->set_cbout(this);
      sb->set_oracle_data(oracle_data);
    }
    for (Integer i = 0; i < psc_dim.dim(); i++) {
      assert(psc_dim(i) > 0);
      use_trace = true;
      modeldim += (psc_dim(i) * (psc_dim(i) + 1) / 2);
      mu_dim += psc_dim(i);
      PSCIPBundleBlock* pb = new PSCIPBundleBlock(psc_dim(i), this);
      pscblock.push_back(pb);
      block.push_back(pb);
      pb->set_cbout(this);
      pb->set_oracle_data(oracle_data);
    }
    if (box_lb.dim() > 0) {
      if ((!scale_box) && (ft != ObjectiveFunction)) {
        if (cb_out())
          get_out() << "**** WARNING: QPConeModelBlock::init(): scale_box is false but the function task is not ObjectionFunction; setting scale_box to true" << std::endl;
        scale_box = true;
      }
      if (modeldim == 0) {
        boxblock = new BoxIPBundleBlock(box_lb, box_ub, scale_box, trace_rhs, this);
        mu_dim += 2 * box_lb.dim() + (scale_box ? 2 : 0);
      } else {
        use_trace |= scale_box;
        boxblock = new BoxIPBundleBlock(box_lb, box_ub, scale_box, -1., this);
        mu_dim += 2 * box_lb.dim() + (scale_box ? 1 : 0);
      }
      boxblock->set_cbout(this);
      modeldim += boxblock->get_vecdim();
      block.push_back(boxblock);
      boxblock->set_oracle_data(oracle_data);
    }
    if ((use_trace) && (ft != ObjectiveFunction))
      mu_dim += 1;

    if (use_trace) {
      //usually the trace vector only needs to be collected once
      trace_vec.init(modeldim, 1, 0.);
      Integer startindex = 0;
      for (unsigned i = 0; i < block.size(); i++) {
        block[i]->add_trace(trace_vec, 1., startindex);
        startindex += block[i]->get_vecdim();
      }
      // std::cout<<" trace_vec="<<trace_vec<<std::endl;
    }

    return 0;
  }

  /// change the right hand side of the trace constraint to b
  int QPConeModelBlock::adjust_trace(Real b) {
    assert(b > 0.); trace_rhs = b; return 0;
  }


  /// evaluate the left hand side of the trace constraint for modelx
  Real QPConeModelBlock::evaluate_trace() const {
    return diff_trace + trace_rhs;
  }

  /// get the right hand side of the trace constraint
  Real QPConeModelBlock::get_trace() {
    return trace_rhs;
  }

  /// get the linear part of modelx (and a guess, which of them are active, in {0.,1.})
  int QPConeModelBlock::get_nncx(Matrix& nncx,
    Matrix* nncx_activity,
    bool cautious) {
    if (nncblock == 0) {
      nncx.init(0, 1, 0.);
      if (nncx_activity)
        nncx_activity->init(0, 1, 0.);
      return 0;
    }
    return nncblock->get_nncx(nncx, nncx_activity, trace_rhs, cautious);
  }

  /// get the SOC part of modelx (and a guess whether the entire cone is active
  int QPConeModelBlock::get_socx(Integer i,
    Matrix& socx,
    Real* socx_activity,
    bool cautious) {
    if (0 > i)
      return 1;
    if (i >= Integer(socblock.size())) {
      socx.init(0, 1, 0.);
      if (socx_activity)
        *socx_activity = 0.;
    }
    assert(socblock[unsigned(i)]);
    return socblock[unsigned(i)]->get_socx(socx, socx_activity, trace_rhs, cautious);
  }

  /// get the PSC part of modelx (and a guess on the rank of the active part)
  int QPConeModelBlock::get_pscx(Integer i,
    Matrix& pscx_eigs,
    Matrix& pscx_vecs,
    Real& pscx_growthrate,
    Matrix& pscx_primalgrowth,
    Matrix& pscx_dualgrowth) {
    if (0 > i)
      return 1;
    if (i >= Integer(pscblock.size())) {
      pscx_eigs.init(0, 1, 0.);
      pscx_vecs.init(0, 1, 0.);
      pscx_growthrate = 0.;
      pscx_primalgrowth.init(0, 1, 0.);
      pscx_dualgrowth.init(0, 1, 0.);
    }
    assert(pscblock[unsigned(i)]);
    return pscblock[unsigned(i)]->get_pscx(pscx_eigs, pscx_vecs,
      pscx_growthrate,
      pscx_primalgrowth, pscx_dualgrowth);
  }

  /// get the linear part of modelx (and a guess, which of them are active, in {0.,1.})
  int QPConeModelBlock::get_boxx(Matrix& boxx,
    Matrix* boxx_activity,
    bool cautious) {
    if (boxblock == 0) {
      boxx.init(0, 1, 0.);
      if (boxx_activity)
        boxx_activity->init(0, 1, 0.);
      return 0;
    }
    return boxblock->get_boxx(boxx, boxx_activity, trace_rhs, cautious);
  }

  /// return the value of the dual variable to the trace consrat == support function value
  Real QPConeModelBlock::tracedual(Real* prec) const {
    if (prec)
      *prec = mu;
    return trace_dual;
  }

  //-----------  QPBlock routines

  /// returns the dimension of the model set (here the same as the bundle size)
  Integer QPConeModelBlock::dim_model() {
    return modeldim;
  }

  /// returns the dimension of the system describing the model set (may contain further constraints)
  Integer QPConeModelBlock::dim_constraints() {
    return use_trace ? 1 : 0;
  }

  ///returns the dual upper bound to the model value (the trace weighted sum of the dual trace variables); it returns 0. if no model is contained
  Real QPConeModelBlock::constraints_cost(void) {
    Real value = trace_dual * trace_rhs;
    if (boxblock)
      value += boxblock->constraints_cost();
    return value;
  }

  //------------ QPConeModelBlockObject routines

  Matrix& QPConeModelBlock::B_times(const Matrix& A,
    Matrix& C,
    Real alpha,
    Real beta,
    int Btrans,
    int Atrans,
    Integer startindex_model,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->B_times(A, C, alpha, beta, Btrans, Atrans, startindex_model,
        globalbundle, startindex_bundle);
      startindex_model += block[i]->get_vecdim();
      startindex_bundle += block[i]->dim_bundle();
    }
    return C;
  }

  /// C=beta*C+alpha*A*B where A and B may be transposed; carry out the model part of this beginning at startindex_model 
  Matrix& QPConeModelBlock::times_B(const Matrix& A,
    Matrix& C,
    Real alpha,
    Real beta,
    int Atrans,
    int Btrans,
    Integer startindex_model,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->times_B(A, C, alpha, beta, Atrans, Btrans, startindex_model,
        globalbundle, startindex_bundle);
      startindex_model += block[i]->get_vecdim();
      startindex_bundle += block[i]->dim_bundle();
    }
    return C;
  }


  ///add the main diagonal block tranpose(projection)*diagvec*projection to bigS starting at startindex
  Symmatrix&
    QPConeModelBlock::add_BDBt(const Matrix& diagvec,
      Symmatrix& bigS,
      bool minus,
      Integer startindex,
      Matrix& Bt,
      Integer startindex_model,
      MinorantBundle& globalbundle,
      Integer startindex_bundle) {
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->add_BDBt(diagvec, bigS, minus, startindex,
        Bt, startindex_model, globalbundle, startindex_bundle);
      startindex_model += block[i]->get_vecdim();
      startindex_bundle += block[i]->dim_bundle();
    }
    return bigS;
  }


  /// get the current matrix for the coupling matrix Bt in the first row of blocks
  Matrix& QPConeModelBlock::get_Bt(Matrix& Bt,
    Integer startindex_model,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->get_Bt(Bt, startindex_model,
        globalbundle, startindex_bundle);
      startindex_model += block[i]->get_vecdim();
      startindex_bundle += block[i]->dim_bundle();
    }
    return Bt;
  }

  /// set the local modelx value in modelx beginning with startindex (initialize it, do not add)
  int QPConeModelBlock::get_modelx(Matrix& inmodelx, Integer startindex) {
    assert(inmodelx.dim() >= startindex + modeldim);
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->vecgetsax(inmodelx, startindex);
      startindex += block[i]->get_vecdim();
    }
    return 0;
  }

  /// set the local modeldx value in modeldx beginning with startindex (initialize it, do not add)
  int QPConeModelBlock::get_modeldx(Matrix& inmodeldx, Integer startindex) {
    assert(inmodeldx.dim() >= startindex + modeldim);
    for (unsigned i = 0; i < block.size(); i++) {
      if (block[i]->get_vecdx(inmodeldx, startindex))
        return 1;
      startindex += block[i]->get_vecdim();
    }
    return 0;
  }

  /// set the local modeldx value in modeldx beginning with startindex (initialize it, do not add)
  int QPConeModelBlock::get_modeldcstr(Matrix& inmodeldcstr, Integer startindex) {
    assert(inmodeldcstr.dim() >= startindex + (use_trace ? 1 : 0));
    inmodeldcstr(startindex) = trace_delta_dual;
    return 0;
  }

  ///returns the aggregate with respect to the current modelx (without constant minorant)
  int QPConeModelBlock::add_modelx_aggregate(Real& val,
    Matrix& vec,
    MinorantBundle& global_bundle,
    Integer startindex_bundle) {
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->add_modelx_aggregate(val, vec, global_bundle, startindex_bundle);
      startindex_bundle += block[i]->dim_bundle();
    }
    return 0;
  }


  /// add  the model violation for the current system solution 
  int QPConeModelBlock::get_sysviol_model(Matrix& insysviol_model,
    Integer startindex_model,
    const Matrix& dy,
    MinorantBundle& global_bundle,
    Integer startindex_bundle) {
    assert(insysviol_model.dim() >= startindex_model + modeldim);
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->get_sysviol_model(insysviol_model, startindex_model,
        dy, trace_delta_dual,
        global_bundle, startindex_bundle);
      startindex_model += block[i]->get_vecdim();
      startindex_bundle += block[i]->dim_bundle();
    }
    return 0;
  }

  /// set the constraint violation for the current system solution 
  int QPConeModelBlock::get_sysviol_constraints(Matrix& constrvec,
    Integer startindex) {
    if ((use_trace) && (modeldim > 0)) {
      assert(constrvec.dim() >= startindex + (use_trace ? 1 : 0));
      Real val = trace_rhs;
      if (ft != ObjectiveFunction)
        val -= trace_slack + trace_delta_slack;
      for (unsigned i = 0; i < block.size(); i++) {
        val -= block[i]->evaluate_trace_x();
        val -= block[i]->evaluate_trace_dx();
      }
      constrvec(startindex) = val;
    }
    return 0;
  }

  void QPConeModelBlock::display_model_values(const Matrix& y,
    MinorantBundle& global_bundle,
    Integer startindex_bundle,
    std::ostream& out) {
    Real val = 0.;
    Matrix vec(y.rowdim(), 1, 0.);
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->add_modelx_aggregate(val, vec, global_bundle, startindex_bundle);
      startindex_bundle += block[i]->dim_bundle();
    }
    out << " modelval=" << val + ip(vec, y) << std::endl;
  }




  // bundlevalues holds the evaluation of the bundle for the current y 
  int QPConeModelBlock::reset_starting_point(const Matrix& y,
    Real mu,
    MinorantBundle& global_bundle,
    Integer startindex_bundle) {
    modelx_changed();

    // //TEST output begin
    // std::cout<<" y="<<y;
    // for (Integer i=startindex_bundle;i<modeldim;i++)
    //   global_bundle[unsigned(i)].display(std::cout);
    // //TEST output end

    //primal side
    Real initval = trace_rhs / mu_dim;
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->center_x(initval);
    }
    if ((!use_trace) || (ft == ObjectiveFunction)) {
      trace_slack = 0;
    } else {
      trace_slack = initval;
    }

    //dual side
    Real add_val = 0.;
    Integer si = startindex_bundle;
    Real tracez = 0.;
    for (unsigned i = 0; i < block.size(); i++) {
      Real av = 0.;
      block[i]->set_bundle_z(y, global_bundle, si, av);
      if (av > add_val)
        add_val = av;
      tracez += block[i]->evaluate_trace_z();
      si += block[i]->dim_bundle();
    }
    if (use_trace) {
      trace_dual = max(add_val + 1., mu / initval - tracez / mu_dim);
      modeldim = 0;
      for (unsigned i = 0; i < block.size(); i++) {
        block[i]->center_z(trace_dual, true);
        block[i]->add_trace_to_diff_model(trace_dual);
        modeldim += block[i]->get_vecdim();
      }

      diff_trace = -trace_slack;
    }
    /*
    diff_model.init(modeldim,1,negbundlevalues.get_store()+startindex);
    diff_model.xpeya(trace_vec,trace_dual);
    */
    return 0;
  }


  ///Euclidean norm of constraint violation of modelx
  Real QPConeModelBlock::primalviol_2normsqr() {
    if (use_trace) {
      Real val = diff_trace + trace_slack;
      return val * val;
    }
    return 0.;
  }

  /// on input viol holds in row i the evalution of bundle minornat i; add the local variables to obtain dual model violation and return in addition to this for the local indices the squared Euclidean norm of this violation  
  Real QPConeModelBlock::dualviol_2normsqr() {
    Real val = 0.;
    for (unsigned i = 0; i < block.size(); i++) {
      val += block[i]->dualviol_2normsqr();
    }
    return val;
  }


  int QPConeModelBlock::get_mu_info(Integer& inmudim,
    Real& tr_xz,
    Real& tr_xdzpdxz,
    Real& tr_dxdz,
    Real& min_xz,
    Real& max_xz) const {
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->get_mu_info(inmudim, tr_xz, tr_xdzpdxz, tr_dxdz, min_xz, max_xz);
    }
    if ((use_trace) && (ft != ObjectiveFunction)) {
      inmudim += 1;
      tr_xz += trace_slack * trace_dual;
      tr_xdzpdxz += (trace_slack * trace_delta_dual) + (trace_delta_slack * trace_dual);
      tr_dxdz += trace_delta_slack * trace_delta_dual;
    }

    return 0;
  }


  int QPConeModelBlock::get_nbh_info(Integer inmudim,
    Real tr_xz,
    Real tr_xdzpdxz,
    Real tr_dxdz,
    Real nbh_ubnd,
    Real& alpha,
    Real& max_nbh,
    Real& nrmsqr_xz,
    Real& nrmsqr_xdzpdxz,
    Real& nrmsqr_dxdz,
    Real& ip_xz_xdzpdxz,
    Real& ip_xz_dxdz,
    Real& ip_dxdz_xdzpdxz) const {
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->get_nbh_info(inmudim, tr_xz, tr_xdzpdxz, tr_dxdz,
        nbh_ubnd, alpha, max_nbh,
        nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
        ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
    }
    if ((use_trace) && (ft != ObjectiveFunction)) {
      const Real mu_xz = tr_xz / inmudim;
      const Real mu_xdzpdxz = tr_xdzpdxz / inmudim;
      const Real mu_dxdz = tr_dxdz / inmudim;
      const Real mu_at_one = mu_xz + mu_xdzpdxz + mu_dxdz;

      NNC_nbh_stepsize(trace_slack, trace_dual, trace_delta_slack, trace_delta_dual,
        mu_xz, mu_xdzpdxz, mu_dxdz, mu_at_one,
        nbh_ubnd, alpha, max_nbh,
        nrmsqr_xz, nrmsqr_xdzpdxz, nrmsqr_dxdz,
        ip_xz_xdzpdxz, ip_xz_dxdz, ip_dxdz_xdzpdxz);
    }

    return 0;
  }


  /// if necessary, reduce alpha to the biggest value so that feasibility is maintained with this step size
  int QPConeModelBlock::linesearch(Real& alpha) const {
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->linesearch(alpha);
    }
    if ((use_trace) && (ft != ObjectiveFunction)) {
      if (trace_delta_slack < -eps_Real)
        alpha = min(alpha, -trace_slack / trace_delta_slack);
      if (trace_delta_dual < -eps_Real)
        alpha = min(alpha, -trace_dual / trace_delta_dual);
    }
    return 0;
  }

  /// compute the step in the model space given the step in the design space
  int QPConeModelBlock::compute_step(const Matrix& ystep,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    modeldx.init(0, 1, 0.);
    //make sure everything is initialized
    //compute sysrhs_model- deltabundle
    //mat_xpey(modeldim,sysrhs_model.get_store(),negdeltabundlevalues.get_store()+startindex);
    Matrix tmpmodelrhs(sysrhs_model);
    //std::cout<<" modrhs="<<tmpmodelrhs;
    B_times(ystep, tmpmodelrhs, -1., 1., 0, 0, 0, globalbundle, startindex_bundle);
    //std::cout<<" modrhsmBdy="<<tmpmodelrhs;
    if (use_trace) {
      if (sysinv_trace.dim() != modeldim) {
        sysinv_trace.init(trace_vec);
        Integer modelind = 0;
        for (unsigned i = 0; i < block.size(); i++) {
          block[i]->apply_xizinv(sysinv_trace, modelind);
          modelind += block[i]->get_vecdim();
        }
      }
      if (ft == ObjectiveFunction) {
        //compute dual delta step
        //std::cout<<" trrhs="<<sysrhs_trace<<" ip=(sysinv_trace,tmpmodelrhs)"<<ip(sysinv_trace,tmpmodelrhs)<<" tr(sysinv_trace)"<<evaluate_trace(sysinv_trace);
        trace_delta_dual = (sysrhs_trace - ip(sysinv_trace, tmpmodelrhs)) / evaluate_trace(sysinv_trace);
        //std::cout<<" trace_delta_dual="<<trace_delta_dual;
        //correct model rhs for this delts step
        tmpmodelrhs.xpeya(trace_vec, trace_delta_dual);
      } else {
        //add trace part to rhs
        tmpmodelrhs.xpeya(trace_vec, trace_dual * sysrhs_trace / trace_slack);
        //add rank one inverse part
        //Real factor=ip(sysinv_trace,sysrhs_model)/(trace_dual/trace_slack+evaluate_trace(sysinv_trace));     
        Real factor = ip(sysinv_trace, tmpmodelrhs) / (trace_slack / trace_dual + evaluate_trace(sysinv_trace));
        tmpmodelrhs.xpeya(trace_vec, -factor);
      }
    }
    //std::cout<<" modrhsmtr="<<tmpmodelrhs;
    tmpmodelrhs *= -1.;
    Integer modelind = 0;
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->set_dx_xizsolverhs(tmpmodelrhs, modelind);
      modelind += block[i]->get_vecdim();
    }
    if ((use_trace) && (ft != ObjectiveFunction)) {
      Real trace_dx = 0.;
      for (unsigned i = 0; i < block.size(); i++)
        trace_dx += block[i]->evaluate_trace_dx();
      trace_delta_dual = trace_dual / trace_slack * (trace_dx + sysrhs_trace);
      trace_delta_slack = -(trace_dx + sysrhs_trace) + complrhs_trace - trace_slack;
      if (std::fabs(diff_trace + trace_slack + trace_delta_slack + trace_dx) > 1e-8 * trace_rhs) {
        if (cb_out()) {
          get_out() << "**** WARNING: QPConeModelBlock::compute_step(...): std::fabs(diff_trace+trace_slack+trace_delta_slack+trace_dx)=" << std::fabs(diff_trace + trace_slack + trace_delta_slack + trace_dx) << "> 1e-8*trace_rhs=" << 1e-8 * trace_rhs << ", numerical difficulties in solving the KKT system might be the cause" << std::endl;
        }
      }
    }

    return 0;
  }


  /// store this computed step locally and compute the missing local dual step information
  int QPConeModelBlock::computed_step(const Matrix& modelxstep,
    Integer startindex_model,
    const Matrix& modelconstrstep,
    Integer startindex_constr) {
    modeldx.init(0, 1, 0.);
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->set_dx(modelxstep, startindex_model);
      startindex_model += block[i]->get_vecdim();
    }
    if (use_trace) {
      trace_delta_dual = modelconstrstep(startindex_constr);
      if (ft != ObjectiveFunction) {
        trace_delta_slack = -trace_delta_dual * trace_slack / trace_dual - trace_slack + complrhs_trace;
      }
    }

    return 0;
  }


  /// move in the last computed step direction by a step of length alpha
  int QPConeModelBlock::do_step(Real alpha,
    const Matrix& y,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    if (use_trace) {
      trace_dual += trace_delta_dual * alpha;
      if (ft != ObjectiveFunction)
        trace_slack += trace_delta_slack * alpha;
    }

    modeldim = 0;
    diff_trace = -trace_rhs;
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->do_bundle_step(alpha, y, globalbundle, startindex_bundle, trace_dual, trace_rhs);
      if (use_trace)
        diff_trace += block[i]->evaluate_trace_x();
      startindex_bundle += block[i]->dim_bundle();
      modeldim += block[i]->get_vecdim();
    }


    if (modeldim != trace_vec.rowdim()) {
      //   //some variables must have been fixed to zero
      //  assert(modeldim<trace_vec.rowdim());
      trace_vec.init(modeldim, 1, 0.);
      Integer si = 0;
      for (unsigned i = 0; i < block.size(); i++) {
        block[i]->add_trace(trace_vec, 1., si);
        si += block[i]->get_vecdim();
      }
    }

    //make sure everything is reinitialized
    modelx_changed();

    old_mu = mu;
    mu = last_rhs_mu;

    return 0;
  }


  int QPConeModelBlock::add_localrhs(Matrix& globalrhs,
    Real rhsmu,
    Real rhscorr,
    Integer startindex_model,
    Integer startindex_constraints,
    bool append,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    //model right hand side
    sysrhs_model.newsize(modeldim, 1); chk_set_init(sysrhs_model, 1);
    Integer si = 0;
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->set_modelrhs(sysrhs_model, rhsmu, rhscorr, si);
      si += block[i]->get_vecdim();
    }

    //trace right hand side
    if (use_trace) {
      last_rhs_mu = rhsmu;
      sysrhs_trace = diff_trace;
      if (ft == ObjectiveFunction) {
        complrhs_trace = 0.;
      } else {
        complrhs_trace = (rhsmu - rhscorr * trace_delta_slack * trace_delta_dual) / trace_dual;
        sysrhs_trace += complrhs_trace;
      }
    }

    //add the right hand sides to the global rhs in the correct way
    if (append) {
      //no further transformation, just add
      mat_xpey(modeldim, globalrhs.get_store() + startindex_model, sysrhs_model.get_store());
      if (use_trace)
        globalrhs(startindex_constraints) += sysrhs_trace;

    } else {
      //add with Schur-complement on bundle to the first block
      Real trsysitr = 0.; //trace of sysinverse*trace
      Matrix tmpvec(sysrhs_model);
      if (use_trace) {
        //make sure everything is initialized by add_BtinvsysB
        assert(sysinv_trace.dim() == trace_vec.dim());
        assert(sys_trace > 0.);
        trsysitr = evaluate_trace(sysinv_trace);
        if (ft != ObjectiveFunction) {
          tmpvec.xpeya(trace_vec, sysrhs_trace / sys_trace);
        }
        //apply first step of "low rank" inverse
        tmpvec.xpeya(trace_vec, ip(sysinv_trace, tmpvec) / (sys_trace + trsysitr));
        //apply as second step the block inverses
      }
      int modelind = 0;
      for (unsigned i = 0; i < block.size(); i++) {
        block[i]->apply_xizinv(tmpvec, modelind);
        modelind += block[i]->get_vecdim();
      }
      // treat the ObjectiveFunction case
      if ((use_trace) && (ft == ObjectiveFunction)) {
        tmpvec.xpeya(sysinv_trace, sysrhs_trace / trsysitr);
      }
      // now aggregate the bundle with these coefficients into the rhs
      /*
      Real dummy;
      for (Integer i=0;i<modeldim;i++){
  globalbundle[unsigned(i+startindex_bundle)].get_minorant(dummy,globalrhs,0,tmpvec(i),true);
      }
      */
      B_times(tmpvec, globalrhs, 1., 1., 1, 0, 0, globalbundle, startindex_bundle);
    }

    return 0;
  }


  ///add the "scaled" minorant outer products to globalsys, where the correct minroants start at the given index
  int QPConeModelBlock::add_BtinvsysB(Symmatrix& globalsys,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    Matrix tmpmat;
    int bundle_ind = startindex_bundle;

    //modelpart
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->add_BtinvsysB(globalsys, globalbundle, bundle_ind);
      bundle_ind += block[i]->dim_bundle();
    }

    //constraintpart
    if (use_trace) {
      if (sys_trace < 0.) {
        sys_trace = 0.;
        if (ft != ObjectiveFunction)
          sys_trace = trace_slack / trace_dual;
      }
      sysinv_trace.init(trace_vec);
      Integer modelind = 0;
      for (unsigned i = 0; i < block.size(); i++) {
        block[i]->apply_xizinv(sysinv_trace, modelind); //need to call this here after the model xiz were computed
        modelind += block[i]->get_vecdim();
      }

      Real factor = 1. / std::sqrt(sys_trace + evaluate_trace(sysinv_trace));
      tmpmat.newsize(globalsys.rowdim(), 1); chk_set_init(tmpmat, 1);
      /*
  for(Integer i=0;i<modeldim;i++){
      Real dummy;
      globalbundle[unsigned(startindex_bundle+i)].get_minorant(dummy,tmpmat,1,sysinv_trace(i)*factor,true);
      }
      */
      B_times(sysinv_trace, tmpmat, factor, 0., 1, 0, 0, globalbundle, startindex_bundle);
      rankadd(tmpmat, globalsys, -1.);
    }

    return 0;
  }

  int QPConeModelBlock::solve_constrsys(const CH_Matrix_Classes::Symmatrix& ABchol,
    const CH_Matrix_Classes::Matrix& LinvABrhs,
    CH_Matrix_Classes::Matrix& LinvABsol,
    CH_Matrix_Classes::Integer startindex_model,
    CH_Matrix_Classes::Matrix& Crhs_and_sol,
    CH_Matrix_Classes::Integer startindex_constraints) {
    int status = 0;
    if (use_trace) {
      Matrix tmpvec(ABchol.rowdim(), 1, 0.);
      mat_xey(modeldim, tmpvec.get_store() + startindex_model, trace_vec.get_store());
      status = ABchol.Chol_Lsolve(tmpvec);
      if (status) {
        if (cb_out())
          get_out() << "**** WARNING QPConeModelBlock::solve_constrsys(......): ABchol.Chol_Lsolve failed and returned " << status << std::endl;
      }
      Real sysval = ip(tmpvec, tmpvec);
      if (ft != ObjectiveFunction)
        sysval += trace_slack / trace_dual;
      assert(sysval > 0.);
      Real rhs = Crhs_and_sol(startindex_constraints) - ip(tmpvec, LinvABrhs);
      rhs /= sysval;
      Crhs_and_sol(startindex_constraints) = rhs;
      LinvABsol.xpeya(tmpvec, rhs);
    }
    return status;
  }



  int QPConeModelBlock::add_localsys(Symmatrix& globalsys,
    Integer startindex_model,
    Integer startindex_constraints) {
    //model block
    if (startindex_model >= 0) {
      Integer modelind = startindex_model;
      for (unsigned i = 0; i < block.size(); i++) {
        block[i]->add_xiz(globalsys, modelind, true);
        modelind += block[i]->get_vecdim();
      }
    }

    //constraint block
    if ((use_trace) && (startindex_constraints >= 0)) {
      for (Integer i = 0; i < modeldim; i++)
        globalsys(startindex_model + i, startindex_constraints) = -trace_vec(i);
      if (sys_trace < 0.) {
        sys_trace = 0.;
        if (ft != ObjectiveFunction)
          sys_trace += trace_slack / trace_dual;
      }
      globalsys(startindex_constraints, startindex_constraints) += sys_trace;
    }

    // //TEST output begin
    // std::cout<<" trace_rhs= "<<trace_rhs;
    // std::cout<<" trace_dual= "<<trace_dual;
    // std::cout<<" trace_slack= "<<trace_slack<<std::endl;
    // std::cout<<"trace_vec= "<<trace_vec<<std::endl;
    // //TEST output end;

    return 0;
  }

  int QPConeModelBlock::localsys_mult(const Matrix& in_vec,
    Matrix& out_vec,
    Integer startindex_model,
    Integer startindex_constraints) {
    //model block
    mat_xey(modeldim, out_vec.get_store() + startindex_model, in_vec.get_store() + startindex_model);
    Integer modelind = startindex_model;
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->apply_xiz(out_vec, modelind, true);
      modelind += block[i]->get_vecdim();
    }
    //constraints part
    if (use_trace) {
      if (sys_trace < 0.) {
        sys_trace = 0.;
        if (ft != ObjectiveFunction)
          sys_trace += trace_slack / trace_dual;
      }
      Real cval = in_vec(startindex_constraints);
      mat_xpeya(modeldim, out_vec.get_store() + startindex_model, trace_vec.get_store(), -cval);

      out_vec(startindex_constraints) = -mat_ip(modeldim, in_vec.get_store() + startindex_model, trace_vec.get_store());
      out_vec(startindex_constraints) += sys_trace * cval;
    }

    return 0;
  }

  /** @brief add the diagonal of the Schur complemented blocks belonging to bundle and local constraints (used for diagonal preconditioning)
  */
  int QPConeModelBlock::add_BCSchur_diagonal(CH_Matrix_Classes::Matrix& diagonal,
    MinorantBundle& globalbundle,
    CH_Matrix_Classes::Integer startindex_bundle) {
    Matrix trafotrace(modeldim, 1, 0.);
    Matrix ipBtrvec(0, 0, 0.);
    Real lambda = 0.;
    if (use_trace) {
      if (sys_trace < 0.) {
        sys_trace = 0.;
        if (ft != ObjectiveFunction)
          sys_trace += trace_slack / trace_dual;
      }

      Integer blockstart_trace = 0;
      for (unsigned i = 0; i < block.size(); i++) {
        block[i]->set_xizinvsqrt_trace(trafotrace, blockstart_trace);
        blockstart_trace += block[i]->dim_bundle();
      }
      Real twonorm = norm2(trafotrace);
      lambda = 1.;
      if (sys_trace > 0.) {
        lambda = 1. / (sys_trace + twonorm * twonorm);
      } else {
        lambda = 1. / (twonorm * twonorm);
      }
      ipBtrvec.init(diagonal.rowdim(), 1, 0.);
    } //endif(use_trace)

    Integer blockstart_trace = 0;
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->add_bundle_xizinv_diagonal(diagonal,
        ipBtrvec,
        globalbundle,
        startindex_bundle,
        trafotrace,
        blockstart_trace);
      blockstart_trace += block[i]->dim_bundle();
      startindex_bundle += block[i]->dim_bundle();
    }

    if (lambda > 0.) {
      for (Integer i = 0; i < diagonal.rowdim(); i++) {
        diagonal(i) -= lambda * sqr(ipBtrvec(i));
      }
    }

    return 0;
  }


  /** @brief append to lowrank "large" columns that should serve well for generating a low rank projection of the Schur complemented model part. For each column i the coordinate sigma_guess(i) gives the Diag_inv-norm for this column. The parameter minval asks to ignore columns whose norms are smaller than minval and minval should be raised to the value indicated by logfraction if the fraction of the log-interval between highest value and 1 is bigger than minval.
  */
  int QPConeModelBlock::propose_BCSchur_pcsubspace(Matrix& lowrank,
    Matrix& sigma_guess,
    const Matrix& Diag_inv,
    Real minval,
    Real diaginvval,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    if (use_trace) {
      //make sure, everything is available
      if (sys_trace < 0.) {
        sys_trace = 0.;
        if (ft != ObjectiveFunction)
          sys_trace += trace_slack / trace_dual;
      }
      if (sysinv_trace.dim() != modeldim) {
        sysinv_trace.init(trace_vec);
        Integer modelind = 0;
        for (unsigned i = 0; i < block.size(); i++) {
          block[i]->apply_xizinv(sysinv_trace, modelind);
          modelind += block[i]->get_vecdim();
        }
      }
      if (schur_trace < 0.) {
        schur_trace = sys_trace + ip(sysinv_trace, trace_vec);
        assert(schur_trace > 0.);
      }
    }
    Real trvectsysinvtrvec = schur_trace - sys_trace;

    Integer si_bundle = startindex_bundle;
    Matrix minus_trmult(0, 1, 0.);
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->add_pcsubspace(lowrank,
        sigma_guess,
        Diag_inv, minval, diaginvval,
        minus_trmult,
        schur_trace,
        globalbundle,
        si_bundle);
      assert(lowrank.coldim() == sigma_guess.rowdim());
      si_bundle += block[i]->dim_bundle();
    }

    if ((schur_trace > 0.) && (norm2(minus_trmult) > 0.)) {
      if (Btsysinv_trace.dim() != lowrank.rowdim()) {
        Btsysinv_trace.init(lowrank.rowdim(), 1, 0.);
        B_times(sysinv_trace, Btsysinv_trace, 1., 0., 1, 0, 0, globalbundle, startindex_bundle);
      }

      Real tracefactor = (1. - std::sqrt(sys_trace / schur_trace)) / trvectsysinvtrvec;
      Indexmatrix delcols;
      for (Integer i = 0; i < minus_trmult.rowdim(); i++) {
        Real d = minus_trmult(i);
        if (d == 0.) {
          continue;
        }
        d *= -tracefactor;
        Integer colind = lowrank.coldim() - minus_trmult.rowdim() + i;
        mat_xpeya(Btsysinv_trace.rowdim(), lowrank.get_store() + colind * lowrank.rowdim(), Btsysinv_trace.get_store(), d);
        if (diaginvval > 0) {
          d = std::sqrt(diaginvval * colip(lowrank, colind));
        } else {
          d = std::sqrt(colip(lowrank, colind, &Diag_inv));
        }
        if (d < minval) {
          delcols.concat_below(colind);
        } else {
          sigma_guess(colind) = d;
        }
      }
      if (delcols.dim() > 0) {
        lowrank.delete_cols(delcols);
        sigma_guess.delete_rows(delcols);
      }
    }

    return 0;
  }

  /** @brief compute the preconditioning low-rank representation of the Schur complementd blocks belonging to bundle and local constraints by adding a Johnson-Lindenstrauss projection onto the given subspace to glob_lowrank

  */
  int QPConeModelBlock::prepare_BCSchur_JLprecond(Matrix& glob_lowrank,
    Matrix& subspace,
    bool append_subspace,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    Matrix* subsp = &subspace;
    Matrix tmpmat;

    if (!append_subspace) {
      Integer startindex_subsp = startindex_bundle;
      if (use_trace) {
        if (sys_trace < 0.) {
          sys_trace = 0.;
          if (ft != ObjectiveFunction)
            sys_trace += trace_slack / trace_dual;
        }
        tmpmat = subspace.rows(Range(startindex_bundle, startindex_bundle + modeldim - 1));
        subsp = &tmpmat;
        startindex_subsp = 0;

        Integer blockstart_trace = 0;
        Matrix trafotrace(modeldim, 1, 0.);
        for (unsigned i = 0; i < block.size(); i++) {
          block[i]->set_xizinvsqrt_trace(trafotrace, blockstart_trace);
          blockstart_trace += block[i]->dim_bundle();
        }

        //eliminate the constraints block by multiplication with  I-lam*v*v'
        Real twonorm = norm2(trafotrace);
        Real lambda = 1.;
        if (sys_trace > 0.) {
          lambda = 1. - std::sqrt(sys_trace / (sys_trace + twonorm * twonorm));
        }
        Matrix tmpvec;
        genmult(trafotrace, tmpmat, tmpvec, 1. / twonorm, 0., 1, 0);
        genmult(trafotrace, tmpvec, tmpmat, -lambda / twonorm, 1.);
      } //endif(use_trace)

      for (unsigned i = 0; i < block.size(); i++) {
        block[i]->add_bundle_xizinvsqrt_projection(glob_lowrank,
          *subsp,
          startindex_subsp,
          globalbundle,
          startindex_bundle);
        startindex_subsp += block[i]->dim_bundle();
        startindex_bundle += block[i]->dim_bundle();
      }
    } else {
      //--- append subspace
      if (use_trace) {
        subsp = &tmpmat;
        tmpmat.newsize(glob_lowrank.coldim(), modeldim);
        tmpmat.init(glob_lowrank.coldim(), 0, 0.);
      }

      for (unsigned i = 0; i < block.size(); i++) {
        block[i]->add_bundle_xizinvsqrt_projection(glob_lowrank,
          *subsp,
          -1,
          globalbundle,
          startindex_bundle);
        startindex_bundle += block[i]->dim_bundle();
      }

      if (use_trace) {
        if (sys_trace < 0.) {
          sys_trace = 0.;
          if (ft != ObjectiveFunction)
            sys_trace += trace_slack / trace_dual;
        }

        Integer blockstart_trace = 0;
        Matrix trafotrace(modeldim, 1, 0.);
        for (unsigned i = 0; i < block.size(); i++) {
          block[i]->set_xizinvsqrt_trace(trafotrace, blockstart_trace);
          blockstart_trace += block[i]->dim_bundle();
        }
        //eliminate the constraints block by multiplication with  I-lam*v*v'
        Real twonorm = norm2(trafotrace);
        Real lambda = 1.;
        if (sys_trace > 0.) {
          lambda = 1. - std::sqrt(1 - twonorm * twonorm / (sys_trace + twonorm * twonorm));
        }
        Matrix tmpvec;
        genmult(tmpmat, trafotrace, tmpvec, 1. / twonorm, 0.);
        genmult(tmpvec, trafotrace, tmpmat, -lambda / twonorm, 1., 0, 1);
        subspace.concat_right(tmpmat);
      }
    } //end else appedn subspace
    return 0;
  }

  /// add the contributions to glob_rhs of the Schur complemented model block, and return local_rhs of the non complemented constraint block in the rows/columns/diagonal block starting at startindex_constraints
  int QPConeModelBlock::add_Schur_rhs(Matrix& glob_rhs,
    Matrix* local_rhs,
    Real rhsmu,
    Real rhscorr,
    Integer startindex_constraints,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    //trace right hand side
    last_rhs_mu = rhsmu;
    if (use_trace) {
      sysrhs_trace = diff_trace;
      if (ft == ObjectiveFunction) {
        complrhs_trace = 0.;
      } else {
        complrhs_trace = (rhsmu - rhscorr * trace_delta_slack * trace_delta_dual) / trace_dual;
        sysrhs_trace += complrhs_trace;
      }
    }

    //collect/add Schur complement of model block
    sysrhs_model.newsize(modeldim, 1); chk_set_init(sysrhs_model, 1);
    Matrix tmpvec(modeldim, 1); chk_set_init(tmpvec, 1);
    Integer modelind = 0;
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->set_modelrhs(sysrhs_model, rhsmu, rhscorr, modelind);
      mat_xey(block[i]->get_vecdim(), tmpvec.get_store() + modelind, sysrhs_model.get_store() + modelind);
      block[i]->apply_xizinv(tmpvec, modelind);
      modelind += block[i]->get_vecdim();
    }

    if (use_trace) {
      if (local_rhs) {
        (*local_rhs)(startindex_constraints) = sysrhs_trace - ip(trace_vec, tmpvec);
      } else {
        //compute the Schur complement of the trace block as well
        if (sys_trace < 0.) {
          sys_trace = 0.;
          if (ft != ObjectiveFunction)
            sys_trace += trace_slack / trace_dual;
        }
        if (sysinv_trace.dim() != modeldim) {
          sysinv_trace.init(trace_vec);
          Integer modelind = 0;
          for (unsigned i = 0; i < block.size(); i++) {
            block[i]->apply_xizinv(sysinv_trace, modelind);
            modelind += block[i]->get_vecdim();
          }
        }
        if (schur_trace < 0.) {
          schur_trace = sys_trace + ip(sysinv_trace, trace_vec);
          assert(schur_trace > 0.);
        }
        Real coeff = (sysrhs_trace - ip(trace_vec, tmpvec)) / schur_trace;
        tmpvec.xpeya(sysinv_trace, coeff); //added below
      }
    }
    B_times(tmpvec, glob_rhs, 1., 1., 1, 0, 0, globalbundle, startindex_bundle);

    return 0;
  }

  /// multiply in_vec with the local contribution to the global main block and add it to out_vec; the other local multiplications are carried out externally with the information provided in prepare_Schur_precond and are not done here.
  int QPConeModelBlock::add_Schur_mult(const Matrix& in_vec,
    Matrix& out_vec,
    const Matrix* in_cvec,
    Matrix* out_cvec,
    Integer startindex_constraints,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    Real cval = 0.;
    Real* out_cval = 0;
    if (use_trace) {
      if (in_cvec) {
        assert(out_cvec);
        cval = (*in_cvec)(startindex_constraints);
        out_cval = &((*out_cvec)(startindex_constraints));
        if (ft != ObjectiveFunction)
          (*out_cval) += cval * trace_slack / trace_dual;
      } else {
        //make sure, everything is available
        if (Btsysinv_trace.dim() != in_vec.dim()) {
          if (sysinv_trace.dim() != modeldim) {
            sysinv_trace.init(trace_vec);
            Integer modelind = 0;
            for (unsigned i = 0; i < block.size(); i++) {
              block[i]->apply_xizinv(sysinv_trace, modelind);
              modelind += block[i]->get_vecdim();
            }
          }
          if (sys_trace < 0.) {
            sys_trace = 0.;
            if (ft != ObjectiveFunction)
              sys_trace += trace_slack / trace_dual;
          }
          if (schur_trace < 0.) {
            schur_trace = sys_trace + ip(sysinv_trace, trace_vec);
            assert(schur_trace > 0.);
          }
          Btsysinv_trace.init(in_vec.dim(), 1, 0.);
          B_times(sysinv_trace, Btsysinv_trace, 1., 0., 1, 0, 0, globalbundle, startindex_bundle);
        }
        cval = ip(Btsysinv_trace, in_vec) / schur_trace;
        //std::cout<<" cval="<<cval<<" schur_trace="<<schur_trace<<std::endl;
      }
    }


    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->add_BtinvsysB_times(in_vec, out_vec,
        cval, out_cval,
        globalbundle, startindex_bundle);
      startindex_bundle += block[i]->dim_bundle();
    }
    return 0;
  }

  int QPConeModelBlock::computed_Schur_step(const Matrix& step_y,
    const Matrix& local_step,
    Integer startindex_constr,
    MinorantBundle& globalbundle,
    Integer startindex_bundle) {
    if (use_trace) {
      trace_delta_dual = local_step(startindex_constr);
      if (ft != ObjectiveFunction) {
        trace_delta_slack = -trace_delta_dual * trace_slack / trace_dual - trace_slack + complrhs_trace;
      }
    } else {
      trace_delta_dual = 0;
      trace_delta_slack = 0.;
    }
    for (unsigned i = 0; i < block.size(); i++) {
      block[i]->set_dx_xizsolvestep(step_y, trace_delta_dual, globalbundle, startindex_bundle);
      startindex_bundle += block[i]->dim_bundle();
    }
    return 0;
  }



}
