#include "Pattern4GPUModule.h"

using namespace Arcane;

/*---------------------------------------------------------------------------*/
/* Message pour la maille en erreur                                          */
/*---------------------------------------------------------------------------*/
void Pattern4GPUModule::
_printError(VariableCellReal cell_arr, const Cell& cell, Real threshold_error) 
{
  std::ostringstream err_msg;

  err_msg 
    << "Cause error:" << std::endl
    << cell_arr.name() << "[" << cell.localId() << "]= " << cell_arr[cell] 
    << " <= " << threshold_error << std::endl;

  throw FatalErrorException(A_FUNCINFO, err_msg.str());
}

/*---------------------------------------------------------------------------*/
/* Implémentation d'origine de computeAndPrintError()                        */
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
_computeAndPrintError_Vori() 
{
  PROF_ACC_BEGIN(__FUNCTION__);

  Real threshold_error = options()->getThresholdError();

  ENUMERATE_CELL(cell_i, allCells()) {
    m_cell_arr1[cell_i] -= m_cell_arr2[cell_i];
    if (m_cell_arr1[cell_i] <= threshold_error) {
      _printError(m_cell_arr1, *cell_i, threshold_error);
    }
  }
  
  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/* Implémentation API GPU Arcane version 0 (sans gestion d'erreur)           */
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
_computeAndPrintError_Varcgpu_v0() 
{
  PROF_ACC_BEGIN(__FUNCTION__);
  
  auto queue = m_acc_env->newQueue();
  auto command = makeCommand(queue);

  auto in_cell_arr2 = ax::viewIn(command, m_cell_arr2);
  auto inout_cell_arr1 = ax::viewInOut(command, m_cell_arr1);

  command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
    inout_cell_arr1[cid] = inout_cell_arr1[cid] - in_cell_arr2[cid];
  };
  
  PROF_ACC_END;
}


/*---------------------------------------------------------------------------*/
/* Implémentation API GPU Arcane version 1                                   */
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
_computeAndPrintError_Varcgpu_v1() 
{
  PROF_ACC_BEGIN(__FUNCTION__);
  
  auto queue = m_acc_env->newQueue();
  auto command = makeCommand(queue);

  auto in_cell_arr2 = ax::viewIn(command, m_cell_arr2);
  auto inout_cell_arr1 = ax::viewInOut(command, m_cell_arr1);

  Real threshold_error = options()->getThresholdError();

  // To get the id of a cell on error
  ax::ReducerMin<Int32> min_cid_on_error(command);
  min_cid_on_error.setValue(std::numeric_limits<Int32>::max());

  command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
    inout_cell_arr1[cid] = inout_cell_arr1[cid] - in_cell_arr2[cid];

    if (inout_cell_arr1[cid] <= threshold_error) {
      // Error: we would like to have the smallest cell id on error
      min_cid_on_error.min(cid.asInt32());
    }
  };
  // The smallest cell id on error if it exists (std::numeric_limits<Int32>::max() otherwise)
  Int32 red_cid_on_error = min_cid_on_error.reduce();

  if (red_cid_on_error < std::numeric_limits<Int32>::max()) {
    // If here, it means it exists a cell on error, we get the smallest one
    Cell cell_on_error(mesh()->itemsInternal(IK_Cell).data(), red_cid_on_error);
    _printError(m_cell_arr1, cell_on_error, threshold_error);
  }
  
  PROF_ACC_END;
}

void Pattern4GPUModule::
_computeAndPrintError_Varcgpu_v2() 
{
  PROF_ACC_BEGIN(__FUNCTION__);
  
  auto queue = m_acc_env->newQueue();
  auto command = makeCommand(queue);

  int*out_min_cid_on_err = min_cid_on_err.data();
  *out_min_cid_on_err = std::numeric_limits<Int32>::max();

  auto in_cell_arr2 = ax::viewIn(command, m_cell_arr2);
  auto inout_cell_arr1 = ax::viewInOut(command, m_cell_arr1);
  Real threshold_error = options()->getThresholdError();

  // To get the id of a cell on error

  command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
    inout_cell_arr1[cid] = inout_cell_arr1[cid] - in_cell_arr2[cid];

    if (inout_cell_arr1[cid] <= threshold_error) {
      // Error: we would like to have the smallest cell id on error
#ifdef ARCCORE_DEVICE_CODE
      atomicMin(out_min_cid_on_err,cid.asInt32());
#else
     *out_min_cid_on_err = cid.asInt32();
#endif
    }
  };
  if (*out_min_cid_on_err < std::numeric_limits<Int32>::max()) {
    // If here, it means it exists a cell on error, we get the smallest one
    Cell cell_on_error(mesh()->itemsInternal(IK_Cell).data(), min_cid_on_err[0]);
    _printError(m_cell_arr1, cell_on_error, threshold_error);
  }
  
  PROF_ACC_END;
}
void Pattern4GPUModule::
_computeAndPrintError_Varcgpu_v3() 
{
  PROF_ACC_BEGIN(__FUNCTION__);
  
  auto queue = m_acc_env->newQueue();
  auto command = makeCommand(queue);

  int*out_min_cid_on_err = min_cid_on_err.data();
  *out_min_cid_on_err = std::numeric_limits<Int32>::max();

  auto in_cell_arr2 = ax::viewIn(command, m_cell_arr2);
  auto inout_cell_arr1 = ax::viewInOut(command, m_cell_arr1);
  Real threshold_error = options()->getThresholdError();

  // To get the id of a cell on error

  command << RUNCOMMAND_ENUMERATE(Cell, cid, allCells()) {
    inout_cell_arr1[cid] = inout_cell_arr1[cid] - in_cell_arr2[cid];

    if (inout_cell_arr1[cid] <= threshold_error) {
      // Error: we would like to have the smallest cell id on error
     *out_min_cid_on_err = cid.asInt32();
    }
  };
  if (*out_min_cid_on_err < std::numeric_limits<Int32>::max()) {
    // If here, it means it exists a cell on error, we get the smallest one
    Cell cell_on_error(mesh()->itemsInternal(IK_Cell).data(), min_cid_on_err[0]);
    _printError(m_cell_arr1, cell_on_error, threshold_error);
  }
  
  PROF_ACC_END;
}
/*---------------------------------------------------------------------------*/
/* Le point d'entrée computeAndPrintError() avec le choix des implémentations */
/*---------------------------------------------------------------------------*/

void Pattern4GPUModule::
computeAndPrintError() {

  PROF_ACC_BEGIN(__FUNCTION__);

  switch (options()->getComputeAndPrintErrorVersion()) {
    case CPEV_ori: _computeAndPrintError_Vori(); break;
    case CPEV_arcgpu_v0: _computeAndPrintError_Varcgpu_v0(); break;
    case CPEV_arcgpu_v1: _computeAndPrintError_Varcgpu_v1(); break;
    case CPEV_arcgpu_v2: _computeAndPrintError_Varcgpu_v2(); break;
    case CPEV_arcgpu_v3: _computeAndPrintError_Varcgpu_v3(); break;
    default: break;
  };

  PROF_ACC_END;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

