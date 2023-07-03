#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*   Type....: Makefile                                                      *
#*   File....: makefile                                                      *
#*   Name....: DFN Pre makefile                                              *
#*   Author..: Thorsten Koch, modified by Christoph Helmberg                 *
#*   Copyright by Authors, All rights reserved                               *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

ARCH		:=	$(shell uname -m | sed -e s/sun../sparc/)
OSTYPE		:=	$(shell uname -s | tr A-Z a-z)
#CXX		=	clang++
#CC		=	clang
CXX		=	g++
CC		=	gcc

#MODE		=       DEBU
MODE		=       OPTI

CONICBUNDLE	=	.
CPPFLAGS	=	-I$(CONICBUNDLE)/include -I$(CONICBUNDLE)/CBsources \
			-I$(CONICBUNDLE)/Matrix -I$(CONICBUNDLE)/Tools -I$(CONICBUNDLE)/cppinterface

CBLIBOBJECT	=	memarray.o CBout.o \
                        MatrixCBSolver.o CBSolver.o \
			CB_CSolver.o CFunction.o cb_cppinterface.o \
                        BundleSolver.o BundleModel.o \
                        BundleWeight.o BundleHKWeight.o BundleRQBWeight.o \
			BundleTerminator.o \
                        Groundset.o GroundsetModification.o \
                        UnconstrainedGroundset.o \
                        LPGroundset.o LPGroundsetModification.o \
                        Minorant.o MinorantUseData.o MinorantPointer.o \
                        BundleData.o SumBlockModel.o SumModel.o \
			SumModelParameters.o \
                        SumBundle.o SumBundleHandler.o \
			SumBundleParametersObject.o SumBundleParameters.o \
                        AffineFunctionTransformation.o AFTModification.o \
                        AFTData.o AFTModel.o \
			ConeModel.o NNCModel.o NNCData.o\
                        NNCBoxSupportFunction.o NNCBoxSupportModification.o \
                        NNCModelParameters.o \
			BoxData.o BoxModel.o BoxModelParameters.o BoxOracle.o \
                        SOCData.o SOCModel.o SOCModelParameters.o \
                        SOCSupportFunction.o SOCSupportModification.o \
                        PSCOracle.o PSCData.o PSCModel.o \
			PSCModelParameters.o PSCVariableMetricSelection.o \
                        PSCAffineFunction.o PSCAffineModification.o \
                        PSCPrimal.o LanczMaxEig.o \
                        Coeffmat.o SparseCoeffmatMatrix.o Bigmatrix.o \
                        ModificationBase.o Modification.o \
                        ModificationTreeData.o \
                        VariableMetric.o VariableMetricSVDSelection.o \
			BundleProxObject.o BundleIdProx.o \
                        BundleDiagonalTrustRegionProx.o \
                        BundleDLRTrustRegionProx.o \
                        BundleLowRankTrustRegionProx.o \
                        BundleDenseTrustRegionProx.o \
			UQPModelBlockObject.o UQPModelBlock.o \
			UQPSumModelBlock.o UQPConeModelBlock.o UQPSolver.o \
			QPModelDataObject.o QPModelBlockObject.o \
			QPSolverObject.o QPSolver.o \
			QPModelBlock.o QPSumModelBlock.o QPConeModelBlock.o \
			InteriorPointBlock.o InteriorPointBundleBlock.o \
			NNCIPBlock.o SOCIPBlock.o PSCIPBlock.o \
			NNCIPBundleBlock.o SOCIPBundleBlock.o PSCIPBundleBlock.o \
			BoxIPBundleBlock.o SOCIPProxBlock.o \
			QPKKTSolverObject.o QPDirectKKTSolver.o \
			QPSolverParameters.o QPSolverBasicStructures.o \
			QPKKTPrecondObject.o QPIterativeKKTSolver.o \
			QPKKTSubspaceHPrecond.o QPIterativeKKTHASolver.o \
			QPIterativeKKTHAeqSolver.o QPKKTSolverComparison.o \
                        indexmat.o matrix.o symmat.o  eigval.o ldl.o chol.o aasen.o \
                        qr.o trisolve.o nnls.o sparssym.o sparsmat.o lanczpol.o \
			IterativeSystemObject.o psqmr.o pcg.o minres.o

CTESTOBJECT	=	c_main.o

CXXTESTOBJECT	=	cxx_main.o

MATTESTOBJECT	=	mat_main.o

MCTOBJECT	=	mc_triangle.o

TARGET		=	lib/libcb.a  t_c t_cxx t_mat mc_triangle

#-----------------------------------------------------------------------------

GCCWARN		=	-W -Wall -pedantic -Wcast-qual -Wwrite-strings \
			-Wnon-virtual-dtor -Wcast-align -Wconversion \
			-Wno-char-subscripts -Wpointer-arith -Wundef \
			-Wno-misleading-indentation -Wno-deprecated-declarations

#--- linux.x86_64.g++ settings ---------------------------------------------------
DEBU.linux.x86_64.g++ =  -g -fPIC
OPTI.linux.x86_64.g++ =  -fPIC -DNDEBUG -O3 -march=native -funroll-loops
WARN.linux.x86_64.g++ =	$(GCCWARN)
DEPD.linux.x86_64.g++ =	-MM
LINK.linux.x86_64.g++ =	-lm
AR.linux.x86_64.g++   =	ar
ARFLAGS.linux.x86_64.g++ =	cr
RANLIB.linux.x86_64.g++ =	ranlib
OPTI.linux.x86_64.gcc =	-DNDEBUG -O3
WARN.linux.x86_64.gcc =	$(GCCWARN)
DEBU.linux.x86_64.gcc = 	-g
DEPD.linux.x86_64.gcc =	-MM
LINK.linux.x86_64.gcc =	-lm
#--- linux.x86_64.clang++ settings ---------------------------------------------------
DEBU.linux.x86_64.clang++ =  -g -std=c++11
OPTI.linux.x86_64.clang++ =  -DNDEBUG  -O3 -march=native -funroll-loops
WARN.linux.x86_64.clang++ =	$(GCCWARN)
DEPD.linux.x86_64.clang++ =	-MM
LINK.linux.x86_64.clang++ =	-lm
AR.linux.x86_64.clang++   =	ar
ARFLAGS.linux.x86_64.clang++ =	cr
RANLIB.linux.x86_64.clang++ =	ranlib
OPTI.linux.x86_64.clang =	-DNDEBUG -O3
WARN.linux.x86_64.clang =	$(GCCWARN)
DEBU.linux.x86_64.clang = 	-g
DEPD.linux.x86_64.clang =	-MM
LINK.linux.x86_64.clang =	-lm
#-----------------------------------------------------------------------------

CXXFLAGS        =       $($(MODE).$(OSTYPE).$(ARCH).$(CXX))\
			$(WARN.$(OSTYPE).$(ARCH).$(CXX))

CCFLAGS		= 	$($(MODE).$(OSTYPE).$(ARCH).$(CC))\
			$(WARN.$(OSTYPE).$(ARCH).$(CC))

LDFLAGS		=	$(LINK.$(OSTYPE).$(ARCH).$(CXX))
DFLAGS		=	$(DEPD.$(OSTYPE).$(ARCH).$(CXX))

AR		=       $(AR.$(OSTYPE).$(ARCH).$(CXX))
ARFLAGS		=	$(ARFLAGS.$(OSTYPE).$(ARCH).$(CXX))
RANLIB		=	$(RANLIB.$(OSTYPE).$(ARCH).$(CXX))

OBJDIR		=	$(MODE).$(OSTYPE).$(ARCH).$(CXX)
OBJCTEST	=	$(addprefix $(OBJDIR)/,$(CTESTOBJECT))
OBJCXXTEST	=	$(addprefix $(OBJDIR)/,$(CXXTESTOBJECT))
OBJMATTEST	=	$(addprefix $(OBJDIR)/,$(MATTESTOBJECT))
OBJMCT		=	$(addprefix $(OBJDIR)/,$(MCTOBJECT))
OBJCBLIB	=	$(addprefix $(OBJDIR)/,$(CBLIBOBJECT))

VPATH	        =       . $(CONICBUNDLE)/Matrix $(CONICBUNDLE)/CBsources $(CONICBUNDLE)/CBtestsources $(CONICBUNDLE)/cppinterface

all:		$(TARGET)

t_c:		$(OBJCTEST) lib/libcb.a
		$(CXX) $(CXXFLAGS) $(OBJCTEST) -Llib -lcb $(LDFLAGS) -o $@

t_cxx:		$(OBJCXXTEST) lib/libcb.a
		$(CXX) $(CXXFLAGS) $(OBJCXXTEST) -Llib -lcb $(LDFLAGS)  -o $@

t_mat:		$(OBJMATTEST) lib/libcb.a
		$(CXX) $(CXXFLAGS) $(OBJMATTEST) -Llib -lcb $(LDFLAGS)  -o $@

mc_triangle:	$(OBJMCT) lib/libcb.a
		$(CXX) $(CXXFLAGS) $(OBJMCT) -Llib -lcb $(LDFLAGS)  -o $@

lib/libcb.a:   	include/CBconfig.hxx $(OBJCBLIB)
		@if [ ! -d lib ]; then mkdir lib; fi
	        $(AR) $(ARFLAGS) lib/libcb.a $(OBJCBLIB)
		$(RANLIB) lib/libcb.a
		$(CXX) -shared -o ../bin/ConicBundle.so $(OBJCBLIB)

clean:
		-rm -rf OPTI.* DEBU.* $(TARGET)

$(OBJDIR)/%.o:	%.cxx
		@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.o:	%.c
		@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
		$(CC) $(CPPFLAGS) $(CCFLAGS) -c $< -o $@

$(OBJDIR)/%.d: %.c
		@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
		$(SHELL) -ec '$(CC) $(DFLAGS) $(CPPFLAGS) $< \
		| sed '\''s/\($*\)\.o[ :]*/$(OBJDIR:/=\/)\/\1.o $(OBJDIR:/=\/)\/\1.d : /g'\'' >$@; \
		[ -s $@ ] || rm -f $@'

$(OBJDIR)/%.d: %.cxx
		@if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi
		$(SHELL) -ec '$(CXX) $(DFLAGS) $(CPPFLAGS) $< \
		| sed '\''s/\($*\)\.o[ :]*/$(OBJDIR:/=\/)\/\1.o $(OBJDIR:/=\/)\/\1.d : /g'\'' >$@; \
		[ -s $@ ] || rm -f $@'

include/CBconfig.hxx: Makefile
ifeq ($(MODE),OPTI)
		rm -f include/CBconfig.hxx
		echo "#ifndef __CBCONFIG_HXX__\n#define __CBCONFIG_HXX__\n\n#define CONICBUNDLE_DEBUG 0\n\n#endif\n" > include/CBconfig.hxx
#		echo -e "#ifndef __CBCONFIG_HXX__\n#define __CBCONFIG_HXX__\n\n#define CONICBUNDLE_DEBUG 0\n\n#endif\n" > include/CBconfig.hxx
else
		rm -f include/CBconfig.hxx
		echo "#ifndef __CBCONFIG_HXX__\n#define __CBCONFIG_HXX__\n\n#define CONICBUNDLE_DEBUG 1\n\n#endif\n" > include/CBconfig.hxx
#		echo -e "#ifndef __CBCONFIG_HXX__\n#define __CBCONFIG_HXX__\n\n#define CONICBUNDLE_DEBUG 1\n\n#endif\n" > include/CBconfig.hxx
endif


#--- for some compilers it is not helpful to generate dependencies
#    automatically, therefore the dependencies are given explicitly

include depend

#include		$(OBJSCB:.o=.d)
#include		$(OBJCBLIB:.o=.d)
