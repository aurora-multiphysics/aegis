#ifndef dynamicBatching__
#define dynamicBatching__

#include <iostream>
#include <vector>
#include "mpi.h"
#include "AegisBase.h"
#include "Inputs.h"

/**
 * Classes for handling MPI dynamic batching
 * Includes three classes - DynamicBatchingBase, Scheduler and Worker 
 * The {Scheduler} sends {Batches} to the {Worker} which contains n {Tasks}
*/

class DynamicBatchingBase : public AegisBase
{
  public:
  DynamicBatchingBase(JsonHandler & aegisParams);

  protected:
  // parameters from config file
  int _batchSize = 16;
  bool _profiling = false;
  bool _debug = false;

  // common member variables used within Scheduler/Worker methods
  MPI_Status _status;
  MPI_Request _request;  
  const int _batchDataTag = 100; // MPI tag for send/recvs for the actual data being sent 
  const int _batchIndexTag = 101; // MPI tag for send/recvs for the index of the work to be carried out
  const int _noMoreWork = -1; // flag to call no more work

  private: 
  void read_params(JsonHandler & aegisParams);
};

class Scheduler : public DynamicBatchingBase
{
  public: 

  explicit Scheduler(JsonHandler & aegisParams, std::vector<double> & totalWorkArray) 
    : DynamicBatchingBase { aegisParams } , _totalWorkArray { totalWorkArray } // inherit base constructor
  {}

  void send_intial_batches();
  void loop_through_batches();
  void full_batch(const int &avaialbleProcess);
  void final_batch(const int &avaialbleProcess);
  std::vector<double> return_total_array();

  protected:

  private:
  const int _numberOfWorkers = (nprocs - 1); // total number of avaialble workers for mpirun
  int _inactiveWorkers = 0; // number of currently inactive workers
  int _BatchesComplete = 0; // number of batches completed
  int _individualTasksComplete = 0; // number of individual tasks complete 
  int _taskIndex = 0;
  std::vector<double> _totalWorkArray; // array to return
};

// class Worker : public DynamicBatchingBase
// {
//   public:
//    explicit Worker(const std::shared_ptr<JsonHandler> & aegisParams) 
//     : DynamicBatchingBase { aegisParams } 
//     {}
  


// };

#endif
