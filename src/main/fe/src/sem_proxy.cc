//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "sem_proxy.h"

#include <cartesian_struct_builder.h>
#include <cartesian_unstruct_builder.h>
#include <sem_solver_acoustic.h>
#include <source_and_receiver_utils.h>

#include <cxxopts.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>

using namespace SourceAndReceiverUtils;

SEMproxy::SEMproxy(const SemProxyOptions& opt)
{
  const int order = opt.order;
  nb_elements_[0] = opt.ex;
  nb_elements_[1] = opt.ey;
  nb_elements_[2] = opt.ez;
  nb_nodes_[0] = opt.ex * order + 1;
  nb_nodes_[1] = opt.ey * order + 1;
  nb_nodes_[2] = opt.ez * order + 1;

  const float spongex = opt.boundaries_size;
  const float spongey = opt.boundaries_size;
  const float spongez = opt.boundaries_size;
  const std::array<float, 3> sponge_size = {spongex, spongey, spongez};
  src_coord_[0] = opt.srcx;
  src_coord_[1] = opt.srcy;
  src_coord_[2] = opt.srcz;

  domain_size_[0] = opt.lx;
  domain_size_[1] = opt.ly;
  domain_size_[2] = opt.lz;

  rcv_coord_[0] = opt.rcvx;
  rcv_coord_[1] = opt.rcvy;
  rcv_coord_[2] = opt.rcvz;

  bool isModelOnNodes = opt.isModelOnNodes;
  isElastic_ = opt.isElastic;
  cout << boolalpha;
  bool isElastic = isElastic_;

  const SolverFactory::methodType methodType = getMethod(opt.method);
  const SolverFactory::implemType implemType = getImplem(opt.implem);
  const SolverFactory::meshType meshType = getMesh(opt.mesh);
  const SolverFactory::modelLocationType modelLocation =
      isModelOnNodes ? SolverFactory::modelLocationType::OnNodes
                     : SolverFactory::modelLocationType::OnElements;
  const SolverFactory::physicType physicType = SolverFactory::physicType::Acoustic;

  float lx = domain_size_[0];
  float ly = domain_size_[1];
  float lz = domain_size_[2];
  int ex = nb_elements_[0];
  int ey = nb_elements_[1];
  int ez = nb_elements_[2];

  if (meshType == SolverFactory::Struct)
  {
    switch (order)
    {
      case 1: {
        model::CartesianStructBuilder<float, int, 1> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      case 2: {
        model::CartesianStructBuilder<float, int, 2> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      case 3: {
        model::CartesianStructBuilder<float, int, 3> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      default:
        throw std::runtime_error(
            "Order other than 1 2 3 is not supported (semproxy)");
    }
  }
  else if (meshType == SolverFactory::Unstruct)
  {
    model::CartesianParams<float, int> param(order, ex, ey, ez, lx, ly, lz,
                                             isModelOnNodes);
    model::CartesianUnstructBuilder<float, int> builder(param);
    m_mesh = builder.getModel();
  }
  else
  {
    throw std::runtime_error("Incorrect mesh type (SEMproxy ctor.)");
  }

  // time parameters
  if (opt.autodt)
  {
    float cfl_factor = (order == 2) ? 0.5 : 0.7;
    dt_ = find_cfl_dt(cfl_factor);
  }
  else
  {
    dt_ = opt.dt;
  }
  timemax_ = opt.timemax;
  num_sample_ = timemax_ / dt_;

  m_solver = SolverFactory::createSolver(methodType, implemType, meshType,
                                         modelLocation, physicType, order);
  m_solver->computeFEInit(*m_mesh, sponge_size, opt.surface_sponge,
                          opt.taper_delta);

  initFiniteElem();

  std::cout << "Number of node is " << m_mesh->getNumberOfNodes() << std::endl;
  std::cout << "Number of element is " << m_mesh->getNumberOfElements()
            << std::endl;
  std::cout << "Launching the Method " << opt.method << ", the implementation "
            << opt.implem << " and the mesh is " << opt.mesh << std::endl;
  std::cout << "Model is on " << (isModelOnNodes ? "nodes" : "elements")
            << std::endl;
  std::cout << "Physics type is " << (isElastic ? "elastic" : "acoustic")
            << std::endl;
  std::cout << "Order of approximation will be " << order << std::endl;
  std::cout << "Time step is " << dt_ << "s" << std::endl;
  std::cout << "Simulated time is " << timemax_ << "s" << std::endl;

  is_snapshots_       = opt.saveSnapshots;
  snap_time_interval_ = opt.snapshotInterval;
  snap_folder_        = opt.snapshotFolder;

  std::cout << "[SNAPSHOTS] enabled=" << std::boolalpha << is_snapshots_
            << " interval=" << snap_time_interval_
            << " folder=" << snap_folder_ << std::endl;

  is_sismos_       = opt.saveSismos;
  sismos_input_file_ = opt.sismosInputFile;
  sismos_folder_  = opt.sismosFolder;

  std::cout << "[SISMOS] enabled=" << std::boolalpha << is_sismos_
            << " input=" << sismos_input_file_
            << " folder=" << sismos_folder_ << std::endl;

}

void SEMproxy::run()
{
  time_point<system_clock> startComputeTime, startOutputTime, totalComputeTime,
      totalOutputTime;

  SEMsolverDataAcoustic solverData(i1, i2, myRHSTerm, pnGlobal, rhsElement,
                                   rhsWeights);


  if (is_sismos_) {
      initSismos();
  }

  for (int indexTimeSample = 0; indexTimeSample < num_sample_;
       indexTimeSample++)
  {
    startComputeTime = system_clock::now();
    m_solver->computeOneStep(dt_, indexTimeSample, solverData);
    totalComputeTime += system_clock::now() - startComputeTime;

    startOutputTime = system_clock::now();

    if (indexTimeSample % 50 == 0)
    {
      m_solver->outputSolutionValues(indexTimeSample, i1, rhsElement[0],
                                     pnGlobal, "pnGlobal");
    }

    // Save pressure at receiver
    const int order = m_mesh->getOrder();

    float varnp1 = 0.0;
    for (int i = 0; i < order + 1; i++)
    {
      for (int j = 0; j < order + 1; j++)
      {
        for (int k = 0; k < order + 1; k++)
        {
          int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[0], i, j, k);
          int globalNodeOnElement =
              i + j * (order + 1) + k * (order + 1) * (order + 1);
          varnp1 +=
              pnGlobal(nodeIdx, i2) * rhsWeightsRcv(0, globalNodeOnElement);
        }
      }
    }

    pnAtReceiver(0, indexTimeSample) = varnp1;

    if (is_sismos_) {
        saveSismos(indexTimeSample);
    }

    if (is_snapshots_ && (indexTimeSample % snap_time_interval_ == 0))
    {
      saveSnapshot(indexTimeSample);
    }

    swap(i1, i2);

    auto tmp = solverData.m_i1;
    solverData.m_i1 = solverData.m_i2;
    solverData.m_i2 = tmp;

    totalOutputTime += system_clock::now() - startOutputTime;
  }

  if (is_sismos_) {
      closeSismos();
  }

  float kerneltime_ms = time_point_cast<microseconds>(totalComputeTime)
                            .time_since_epoch()
                            .count();
  float outputtime_ms =
      time_point_cast<microseconds>(totalOutputTime).time_since_epoch().count();

  cout << "------------------------------------------------ " << endl;
  cout << "\n---- Elapsed Kernel Time : " << kerneltime_ms / 1E6 << " seconds."
       << endl;
  cout << "---- Elapsed Output Time : " << outputtime_ms / 1E6 << " seconds."
       << endl;
  cout << "------------------------------------------------ " << endl;


    //4. MESURE DU TEMPS
  saveMetrics(totalComputeTime, totalOutputTime);
}

// Initialize arrays
void SEMproxy::init_arrays()
{
  cout << "Allocate host memory for source and pressure values ..." << endl;

  rhsElement = allocateVector<vectorInt>(myNumberOfRHS, "rhsElement");
  rhsWeights = allocateArray2D<arrayReal>(
      myNumberOfRHS, m_mesh->getNumberOfPointsPerElement(), "RHSWeight");
  myRHSTerm = allocateArray2D<arrayReal>(myNumberOfRHS, num_sample_, "RHSTerm");
  pnGlobal =
      allocateArray2D<arrayReal>(m_mesh->getNumberOfNodes(), 2, "pnGlobal");
  pnAtReceiver = allocateArray2D<arrayReal>(1, num_sample_, "pnAtReceiver");
  // Receiver
  rhsElementRcv = allocateVector<vectorInt>(1, "rhsElementRcv");
  rhsWeightsRcv = allocateArray2D<arrayReal>(
      1, m_mesh->getNumberOfPointsPerElement(), "RHSWeightRcv");
}

// Initialize sources
void SEMproxy::init_source()
{
  arrayReal myRHSLocation = allocateArray2D<arrayReal>(1, 3, "RHSLocation");
  // std::cout << "All source are currently are coded on element 50." <<
  // std::endl;
  std::cout << "All source are currently are coded on middle element."
            << std::endl;
  int ex = nb_elements_[0];
  int ey = nb_elements_[1];
  int ez = nb_elements_[2];

  int lx = domain_size_[0];
  int ly = domain_size_[1];
  int lz = domain_size_[2];

  // Get source element index

  int source_index = floor((src_coord_[0] * ex) / lx) +
                     floor((src_coord_[1] * ey) / ly) * ex +
                     floor((src_coord_[2] * ez) / lz) * ey * ex;

  for (int i = 0; i < 1; i++)
  {
    rhsElement[i] = source_index;
  }

  // Get coordinates of the corners of the sourc element
  float cornerCoords[8][3];
  int I = 0;
  int nodes_corner[2] = {0, m_mesh->getOrder()};
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh->globalNodeIndex(rhsElement[0], i, j, k);
        cornerCoords[I][0] = m_mesh->nodeCoord(nodeIdx, 0);
        cornerCoords[I][2] = m_mesh->nodeCoord(nodeIdx, 2);
        cornerCoords[I][1] = m_mesh->nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

  // initialize source term
  vector<float> sourceTerm =
      myUtils.computeSourceTerm(num_sample_, dt_, f0, sourceOrder);
  for (int j = 0; j < num_sample_; j++)
  {
    myRHSTerm(0, j) = sourceTerm[j];
    if (j % 100 == 0)
      cout << "Sample " << j << "\t: sourceTerm = " << sourceTerm[j] << endl;
  }

  // get element number of source term
  myElementSource = rhsElement[0];
  cout << "Element number for the source location: " << myElementSource << endl
       << endl;

  int order = m_mesh->getOrder();

  switch (order)
  {
    case 1:
      SourceAndReceiverUtils::ComputeRHSWeights<1>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    case 2:
      SourceAndReceiverUtils::ComputeRHSWeights<2>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    case 3:
      SourceAndReceiverUtils::ComputeRHSWeights<3>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    default:
      throw std::runtime_error("Unsupported order: " + std::to_string(order));
  }

  // Receiver computation
  int receiver_index = floor((rcv_coord_[0] * ex) / lx) +
                       floor((rcv_coord_[1] * ey) / ly) * ex +
                       floor((rcv_coord_[2] * ez) / lz) * ey * ex;

  for (int i = 0; i < 1; i++)
  {
    rhsElementRcv[i] = receiver_index;
  }

  // Get coordinates of the corners of the receiver element
  float cornerCoordsRcv[8][3];
  I = 0;
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[0], i, j, k);
        cornerCoordsRcv[I][0] = m_mesh->nodeCoord(nodeIdx, 0);
        cornerCoordsRcv[I][2] = m_mesh->nodeCoord(nodeIdx, 2);
        cornerCoordsRcv[I][1] = m_mesh->nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

  switch (order)
  {
    case 1:
      SourceAndReceiverUtils::ComputeRHSWeights<1>(cornerCoordsRcv, rcv_coord_,
                                                   rhsWeightsRcv);
      break;
    case 2:
      SourceAndReceiverUtils::ComputeRHSWeights<2>(cornerCoordsRcv, rcv_coord_,
                                                   rhsWeightsRcv);
      break;
    case 3:
      SourceAndReceiverUtils::ComputeRHSWeights<3>(cornerCoordsRcv, rcv_coord_,
                                                   rhsWeightsRcv);
      break;
    default:
      throw std::runtime_error("Unsupported order: " + std::to_string(order));
  }
}

SolverFactory::implemType SEMproxy::getImplem(string implemArg)
{
  if (implemArg == "makutu") return SolverFactory::MAKUTU;
  if (implemArg == "shiva") return SolverFactory::SHIVA;

  throw std::invalid_argument(
      "Implentation type does not follow any valid type.");
}

SolverFactory::meshType SEMproxy::getMesh(string meshArg)
{
  if (meshArg == "cartesian") return SolverFactory::Struct;
  if (meshArg == "ucartesian") return SolverFactory::Unstruct;

  std::cout << "Mesh type found is " << meshArg << std::endl;
  throw std::invalid_argument("Mesh type does not follow any valid type.");
}

SolverFactory::methodType SEMproxy::getMethod(string methodArg)
{
  if (methodArg == "sem") return SolverFactory::SEM;
  if (methodArg == "dg") return SolverFactory::DG;

  throw std::invalid_argument("Method type does not follow any valid type.");
}

float SEMproxy::find_cfl_dt(float cfl_factor)
{
  float sqrtDim3 = 1.73;  // to change for 2d
  float min_spacing = m_mesh->getMinSpacing();
  float v_max = m_mesh->getMaxSpeed();

  float dt = cfl_factor * min_spacing / (sqrtDim3 * v_max);

  return dt;
}

void SEMproxy::saveSnapshot(int timestep)
{
  if (!is_snapshots_) return;

  std::filesystem::create_directories(snap_folder_);

  std::ostringstream oss;
  oss << snap_folder_ << "/snapshot_"
      << std::setw(6) << std::setfill('0') << timestep << ".dat";
  std::string filename = oss.str();

  std::ofstream out(filename);
  if (!out)
  {
    std::cerr << "Error: cannot open " << filename << " for writing\n";
    return;
  }

  int nNodes = static_cast<int>(m_mesh->getNumberOfNodes());
  float time = timestep * dt_;  // dt_ est le pas de temps stocké dans SEMproxy

  for (int n = 0; n < nNodes; ++n)
  {
    float x = m_mesh->nodeCoord(n, 0);  // 0 = X
    float y = m_mesh->nodeCoord(n, 1);  // 1 = Y
    float z = m_mesh->nodeCoord(n, 2);  // 2 = Z


    float p = pnGlobal(n, i2);

    out << x << " " << y << " " << z << " " << p << "\n";
  }

  out.close();
}

// Initialisation : Lecture du CSV et création des fichiers
void SEMproxy::initSismos()
{
    if (!is_sismos_) return;

    std::cout << "[SISMOS] Initializing seismograms from: " << sismos_input_file_ << std::endl;

    // 1. Création du dossier de sortie
    std::filesystem::create_directories(sismos_folder_);

    // 2. Ouverture du fichier d'entrée CSV
    std::ifstream infile(sismos_input_file_);
    if (!infile.is_open())
    {
        std::cerr << "Error: Cannot open sismos input file: " << sismos_input_file_ << std::endl;
        return;
    }

    std::string line;
    while (std::getline(infile, line))
    {
        if (line.empty()) continue;

        try 
        {
            int nodeID = std::stoi(line);
            
            // Vérification que le noeud est valide dans le maillage
            if (nodeID < 0 || nodeID >= m_mesh->getNumberOfNodes()) {
                std::cerr << "Warning: Node ID " << nodeID << " is out of bounds. Skipped.\n";
                continue;
            }

            // Récupération des coordonnées pour le nom du fichier
            float x = m_mesh->nodeCoord(nodeID, 0);
            float y = m_mesh->nodeCoord(nodeID, 1);
            float z = m_mesh->nodeCoord(nodeID, 2);

            // Construction du nom de fichier : sismo_X_Y_Z.dat
            std::ostringstream oss;
            oss << sismos_folder_ << "/sismo_x" 
                << std::fixed << std::setprecision(2) << x << "_y" 
                << y << "_z" << z << ".csv";
            
            std::string filename = oss.str();

            // Création et ouverture du fichier de sortie
            auto outfile = std::make_unique<std::ofstream>(filename);
            if (!outfile->is_open()) {
                std::cerr << "Error: Cannot create output file " << filename << "\n";
                continue;
            }

            // Écriture d'un header optionnel (Temps Pression)
            *outfile << "Time,Pressure\n";

            // Stockage de l'ID et du pointeur vers le fichier
            sismos_node_ids_.push_back(nodeID);
            sismos_files_.push_back(std::move(outfile));

        }
        catch (const std::exception& e)
        {
            std::cerr << "Error parsing line '" << line << "': " << e.what() << std::endl;
        }
    }
    
    std::cout << "[SISMOS] " << sismos_node_ids_.size() << " probes initialized." << std::endl;
    infile.close();
}

// Écriture : Appelée à chaque pas de temps
void SEMproxy::saveSismos(int timestep)
{
    if (sismos_node_ids_.empty()) return;

    float time = timestep * dt_;

    // Pour chaque noeud observé
    for (size_t i = 0; i < sismos_node_ids_.size(); ++i)
    {
        int nodeID = sismos_node_ids_[i];
        
        // Récupération de la pression (attention à utiliser i2 comme dans saveSnapshot)
        float pressure = pnGlobal(nodeID, i2); 

        // Écriture dans le fichier correspondant
        // Format: Temps [espace] Pression
        *sismos_files_[i] << time << "," << pressure << "\n";
    }
}

// Nettoyage : Fermeture des fichiers
void SEMproxy::closeSismos()
{
    if (!is_sismos_) return;

    for (auto& file : sismos_files_)
    {
        if (file && file->is_open()) {
            file->close();
        }
    }
    sismos_files_.clear();
    sismos_node_ids_.clear();
    std::cout << "[SISMOS] Files closed." << std::endl;
}






// TIME EXECUTION FOR AD HOC
std::string getCurrentDateTime() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);

    std::tm tm{};
#ifdef _WIN32
    localtime_s(&tm, &t);   // safe windows
#else
    localtime_r(&t, &tm);   // safe Linux
#endif

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d_%H%M%S");  // compact format: 20251121_081230
    return oss.str();
}
void SEMproxy::saveMetrics(std::chrono::system_clock::time_point compute_tp,
                           std::chrono::system_clock::time_point output_tp)
{
  
    namespace fs = std::filesystem;

    // Create folder build/metrics if it does not exist
    fs::path metricsDir = "snapshot_metrics";
    if (!fs::exists(metricsDir)) {
        fs::create_directory(metricsDir);
    }

    fs::path metricsFile = metricsDir / "metrics.json";
    //fs::path metricsDir = fs::current_path() / "snapshot-metrics";

    // Duration in seconds
    auto compute_s = std::chrono::duration_cast<std::chrono::microseconds>(compute_tp.time_since_epoch()).count() / 1e6;
    auto output_s  = std::chrono::duration_cast<std::chrono::microseconds>(output_tp.time_since_epoch()).count() / 1e6;

    size_t mesh_nodes    = m_mesh->getNumberOfNodes();
    size_t mesh_elements = m_mesh->getNumberOfElements();


    // We open the file in Read/Write mode
    std::fstream out(metricsFile, std::ios::in | std::ios::out);

    bool first = true;

    if (!fs::exists(metricsFile) || fs::file_size(metricsFile) == 0)
    {
        // New file, we open the object array
        out.close();
        out.open(metricsFile, std::ios::out);
        out << "[\n";
        first = true;
    }
    else
    {
        // File existing: move ] before the close of the array
        out.seekp(-2, std::ios::end); // move pointer
        out.write("", 0);             // nothing is written, just to stablish the position
        out.flush();                  // make sure the pointer is correct
        first = false;
    }



    out << std::fixed << std::setprecision(6);
    out << "  ,{\n";
    out << "    \"date\": \""  << getCurrentDateTime() << "\",\n";
    out << "    \"compute_time\": " << compute_s << ",\n";
    out << "    \"output_time\": " << output_s << ",\n";
    out << "    \"nodes\": " << mesh_nodes << ",\n";
    out << "    \"elements\": " << mesh_elements << ",\n";
    out << "    \"samples\": " << num_sample_ << ",\n";
    out << "    \"snapshot_enabled\": " << is_snapshots_ << ",\n";
    out << "    \"snapshot_interval\": " << snap_time_interval_ << "\n";
    out << "  }\n"; //snapshot_size

    out << "]\n"; // we close the object array
    out.close();
}