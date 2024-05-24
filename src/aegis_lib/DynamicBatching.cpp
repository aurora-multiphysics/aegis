#include "DynamicBatching.h"

DynamicBatchingBase::DynamicBatchingBase(JsonHandler & aegisParams)
{
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  read_params(aegisParams);
}

void
DynamicBatchingBase::read_params(JsonHandler & aegisParams)
{
  if (aegisParams.data().contains("dynamic_batching_params"))
  {
    JsonHandler dynamicBatchingParams(aegisParams.data()["dynmaic_batching_params"]);
    _batchSize = dynamicBatchingParams.get_optional<int>("batch_size").value_or(_batchSize);
    _profiling = dynamicBatchingParams.get_optional<bool>("profling").value_or(_profiling);
    _debug = dynamicBatchingParams.get_optional<bool>("debug").value_or(_debug);
  }
  // std::cout << "MPI Dynamic batching running in debug mode. Data will be written out to files:
  // {DynamicBatchingBase.txt, rank_1.txt ... rank_n.txt} " << std::endl; std::cout << "MPI Dynamic
  // batching running with worker profling enabled. Timings will be printed to stdout..." <<
  // std::endl;
}

// send the initial batches to each process. Perhaps can be rolled into the main loop?
void
Scheduler::send_intial_batches()
{
  int intialWorkerIndexes = 0;
  std::vector<double> workerValues(_batchSize);
  std::cout << "Dynamic batch scheduling with " << _numberOfWorkers << " processes, each handling "
            << _batchSize << " tasks" << std::endl;

  for (int procID = 1; procID < nprocs; procID++)
  {
    MPI_Send(&_individualTasksComplete, 1, MPI_INT, procID, procID, MPI_COMM_WORLD);
    _individualTasksComplete += _batchSize;
  }
}

void
Scheduler::loop_through_batches()
{
  do
  {
    int avaialbleProcess = 0;
    // get available process
    MPI_Recv(&avaialbleProcess, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &_status);

    // if there is still more than a full batch worth of work send entire batch
    if (_individualTasksComplete < _totalWorkArray.size())
    {
      full_batch(avaialbleProcess);
    }
    else
    {
      final_batch(avaialbleProcess);
    }

  } while (_inactiveWorkers < nprocs - 1); // while some workers are still active
}

// Send the work indexes and recieve the full batch size from workers
void
Scheduler::full_batch(const int & avaialbleProcess)
{
  std::vector<double> workerValues(_batchSize);

  MPI_Isend(&_individualTasksComplete, 1, MPI_INT, avaialbleProcess, 1, MPI_COMM_WORLD, &_request);

  MPI_Irecv(workerValues.data(), workerValues.size(), MPI_DOUBLE, MPI_ANY_SOURCE, _batchDataTag,
            MPI_COMM_WORLD, &_request);
  MPI_Irecv(&_taskIndex, 1, MPI_INT, MPI_ANY_SOURCE, _batchIndexTag, MPI_COMM_WORLD, &_request);

  _individualTasksComplete += _batchSize;
  MPI_Wait(&_request, &_status);
  double * totalWorkPtr = _totalWorkArray.data() + _taskIndex;
  memcpy(totalWorkPtr, _totalWorkArray.data(), sizeof(double) * _totalWorkArray.size());
}

// Send message to workers to tell them no more work and recieve final work array back from them
void
Scheduler::final_batch(const int & avaialbleProcess)
{
  // send message to workers to tell them no more work avaialble
  MPI_Send(&_noMoreWork, 1, MPI_INT, avaialbleProcess, 1, MPI_COMM_WORLD);
  _inactiveWorkers++;

  // get size of incoming final worker arrays
  int finalArraySize = 0;
  MPI_Probe(MPI_ANY_SOURCE, _batchDataTag, MPI_COMM_WORLD, &_status);
  MPI_Get_count(&_status, MPI_DOUBLE, &finalArraySize);
  std::vector<double> finalWorkArray(finalArraySize);

  // recieve final worker arrays
  MPI_Recv(finalWorkArray.data(), finalWorkArray.size(), MPI_DOUBLE, _status.MPI_SOURCE,
           _batchDataTag, MPI_COMM_WORLD, &_status);
  MPI_Recv(&_taskIndex, 1, MPI_INT, _status.MPI_SOURCE, _batchIndexTag, MPI_COMM_WORLD, &_status);

  double * totalWorkPtr = _totalWorkArray.data() + _taskIndex;
  memcpy(totalWorkPtr, finalWorkArray.data(), sizeof(double) * finalWorkArray.size());
}

std::vector<double>
Scheduler::return_total_array()
{
  return _totalWorkArray;
}