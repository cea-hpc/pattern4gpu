#ifndef ACC_ENV_DEFAULT_OPTIONS_H
#define ACC_ENV_DEFAULT_OPTIONS_H

enum eDeviceAffinity {
  DA_none = 0,      //! Aucune sélection de device
  DA_cu_world_rank, //! device = world_rank%cuda_device_count
  DA_cu_node_rank,  //! device =  node_rank%cuda_device_count
  DA_cu_hwloc       //! cherche à placer sur le device le plus proche (si hwloc détecté)
};

#endif
