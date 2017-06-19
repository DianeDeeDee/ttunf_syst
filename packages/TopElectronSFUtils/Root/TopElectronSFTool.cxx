#include "TopElectronSFUtils/TopElectronSFTool.h"


const double TopElectronSFTool::ele_ID_etamax = 2.47;

const double TopElectronSFTool::ele_ID_etacrack[2] = {1.37,1.52};

const double TopElectronSFTool::ele_ID_EtMin = 20.;

const double TopElectronSFTool::ele_ID_etabins[18] = {-2.47,-2.37,-2.01,-1.81,-1.37,-1.15,-0.8,-0.6,-0.1,0,0.1,0.6,0.8,1.15,1.52,1.81,2.01,2.37};

const double TopElectronSFTool::ele_ID_ETbins[8] = {20.,25.,30.,35.,40.,45., 50., 60. }; // lower edges of et bins

const double TopElectronSFTool::ele_ID_SFmatrix[8][18] = {  // TightPP + isolation
  {0.918, 0.924, 0.928, 0.958, 0.960, 0.954, 0.932, 0.937, 0.942, 0.936, 0.934, 0.936, 0.951, 0.929, 0.974, 0.949, 0.915, 0.903},
  {0.945, 0.948, 0.964, 0.990, 0.986, 0.986, 0.961, 0.964, 0.977, 0.975, 0.966, 0.963, 0.985, 0.962, 1.000, 0.971, 0.944, 0.949},
  {0.973, 0.964, 0.977, 1.005, 1.004, 1.005, 0.982, 0.984, 0.999, 0.990, 0.980, 0.977, 0.991, 0.967, 1.016, 0.983, 0.962, 0.959},
  {0.964, 0.968, 0.990, 1.016, 1.000, 1.006, 0.990, 0.992, 1.003, 1.006, 0.988, 0.984, 1.000, 0.953, 1.037, 0.985, 0.969, 0.964},
  {0.960, 0.970, 0.979, 1.023, 0.985, 1.003, 0.984, 0.990, 1.003, 1.000, 0.987, 0.981, 1.000, 0.930, 1.041, 0.995, 0.963, 0.962},
  {0.944, 0.957, 0.978, 1.018, 0.987, 0.993, 0.981, 0.983, 0.996, 0.987, 0.984, 0.979, 0.995, 0.913, 1.035, 0.987, 0.954, 0.950},
  {0.938, 0.954, 0.977, 1.018, 0.977, 0.987, 0.976, 0.979, 0.997, 0.994, 0.978, 0.967, 0.983, 0.913, 1.026, 0.976, 0.936, 0.945},
  {0.962, 0.954, 0.973, 1.018, 0.974, 0.990, 0.988, 0.965, 0.993, 0.987, 0.981, 0.968, 0.980, 0.910, 1.020, 0.938, 0.936, 0.940}};
    
const double TopElectronSFTool::ele_ID_errmatrix[8][18] = { // TightPP + isolation
  {0.034, 0.031, 0.032, 0.032, 0.030, 0.030, 0.030, 0.029, 0.036, 0.036, 0.029, 0.030, 0.030, 0.031, 0.033, 0.033, 0.034, 0.036},
  {0.027, 0.023, 0.024, 0.024, 0.021, 0.021, 0.021, 0.020, 0.029, 0.029, 0.020, 0.021, 0.021, 0.022, 0.025, 0.025, 0.027, 0.031},
  {0.027, 0.023, 0.023, 0.023, 0.020, 0.020, 0.020, 0.019, 0.029, 0.029, 0.019, 0.020, 0.020, 0.021, 0.025, 0.025, 0.027, 0.031},
  {0.028, 0.023, 0.024, 0.024, 0.021, 0.020, 0.020, 0.020, 0.029, 0.029, 0.020, 0.020, 0.020, 0.021, 0.025, 0.025, 0.028, 0.031},
  {0.027, 0.023, 0.023, 0.024, 0.020, 0.020, 0.020, 0.019, 0.029, 0.029, 0.019, 0.020, 0.020, 0.021, 0.025, 0.025, 0.027, 0.032},
  {0.029, 0.023, 0.024, 0.024, 0.021, 0.020, 0.020, 0.020, 0.029, 0.029, 0.020, 0.020, 0.020, 0.023, 0.026, 0.025, 0.027, 0.033},
  {0.028, 0.023, 0.024, 0.024, 0.020, 0.020, 0.020, 0.019, 0.029, 0.029, 0.019, 0.020, 0.020, 0.022, 0.025, 0.025, 0.027, 0.032},
  {0.034, 0.025, 0.026, 0.026, 0.022, 0.021, 0.021, 0.020, 0.031, 0.031, 0.020, 0.021, 0.021, 0.025, 0.027, 0.027, 0.029, 0.034}};


// Reconstruction
const double TopElectronSFTool::ele_reco_etamax = 2.47;

const double TopElectronSFTool::ele_reco_etacrack[2] = {1.37,1.52};

const double TopElectronSFTool::ele_reco_etabins[9] = {-2.47,-2.01,-1.37,-0.8,-0.1, 0.1,0.8,1.52,2.01};

const double TopElectronSFTool::ele_reco_SFmatrix[9] = { 1.02, 1.007, 1.002, 0.994, 0.992, 0.993, 1.001, 1.006, 1.023};

const double TopElectronSFTool::ele_reco_errmatrix[9] = {0.007, 0.006, 0.006, 0.011, 0.012, 0.010, 0.006, 0.006, 0.007};


// Trigger
const double TopElectronSFTool::ele_trigger_etamax = 2.47;

const double TopElectronSFTool::ele_trigger_etacrack[2] = {1.37,1.52};

const double TopElectronSFTool::ele_trigger_EtMin= 21.;

const double TopElectronSFTool::ele_trigger_etabins[18] = {-2.47,-2.37,-2.01,-1.81,-1.37,-1.15,-0.8,-0.6,-0.1,0,0.1,0.6,0.8,1.15,1.52,1.81,2.01,2.37};

const double TopElectronSFTool::ele_trigger_ETbins[6] = {21.,23.,25,30.,35.,40.}; // lower edges of et bins

const double TopElectronSFTool::ele_trigger_SFmatrix_e20_medium[6][18] = {
  {0.985, 0.972, 0.977, 0.973, 0.991, 0.985, 0.972, 0.998, 0.915, 0.978, 1.010, 1.007, 1.003, 1.008, 0.978, 0.972, 0.979, 0.989},
  {0.981, 0.968, 0.973, 0.970, 0.987, 0.981, 0.969, 0.995, 0.912, 0.974, 1.006, 1.003, 1.000, 1.005, 0.974, 0.968, 0.975, 0.985},
  {0.988, 0.974, 0.979, 0.976, 0.994, 0.987, 0.975, 1.001, 0.917, 0.980, 1.013, 1.010, 1.006, 1.011, 0.981, 0.975, 0.981, 0.992},
  {0.989, 0.976, 0.981, 0.977, 0.995, 0.989, 0.977, 1.003, 0.919, 0.982, 1.014, 1.011, 1.008, 1.013, 0.982, 0.976, 0.983, 0.993},
  {0.991, 0.977, 0.982, 0.979, 0.996, 0.990, 0.978, 1.004, 0.920, 0.983, 1.015, 1.012, 1.009, 1.014, 0.983, 0.977, 0.984, 0.994},
  {0.991, 0.978, 0.983, 0.979, 0.997, 0.991, 0.979, 1.005, 0.921, 0.984, 1.016, 1.013, 1.010, 1.015, 0.984, 0.978, 0.985, 0.995}};

const double TopElectronSFTool::ele_trigger_SFmatrix_e22_medium[6][18] = {
  {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
  {1.022, 0.971, 0.974, 0.968, 0.987, 0.982, 0.971, 0.994, 0.911, 1.007, 1.040, 1.038, 1.037, 1.035, 0.973, 0.965, 0.978, 1.004},
  {1.026, 0.974, 0.978, 0.972, 0.991, 0.986, 0.975, 0.998, 0.915, 1.011, 1.044, 1.041, 1.041, 1.039, 0.977, 0.968, 0.982, 1.008},
  {1.032, 0.980, 0.983, 0.977, 0.996, 0.991, 0.980, 1.004, 0.920, 1.017, 1.050, 1.047, 1.046, 1.044, 0.982, 0.974, 0.987, 1.013},
  {1.032, 0.980, 0.983, 0.977, 0.996, 0.991, 0.980, 1.004, 0.919, 1.017, 1.050, 1.047, 1.046, 1.044, 0.982, 0.973, 0.987, 1.013},
  {1.032, 0.979, 0.983, 0.977, 0.996, 0.991, 0.980, 1.003, 0.919, 1.017, 1.049, 1.047, 1.046, 1.044, 0.982, 0.973, 0.987, 1.013}};

const double TopElectronSFTool::ele_trigger_SFmatrix_e22vh_medium[6][18] = {
  {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
  {0.934, 0.956, 0.962, 0.966, 0.984, 0.985, 0.975, 0.995, 0.934, 0.967, 0.989, 0.995, 0.986, 0.982, 0.967, 0.955, 0.963, 0.914},
  {0.939, 0.960, 0.966, 0.970, 0.989, 0.989, 0.980, 1.000, 0.938, 0.971, 0.994, 0.999, 0.990, 0.987, 0.971, 0.960, 0.968, 0.918},
  {0.941, 0.962, 0.968, 0.973, 0.991, 0.992, 0.982, 1.002, 0.940, 0.974, 0.996, 1.002, 0.993, 0.989, 0.974, 0.962, 0.970, 0.921},
  {0.943, 0.965, 0.971, 0.975, 0.994, 0.994, 0.985, 1.005, 0.943, 0.976, 0.999, 1.004, 0.995, 0.992, 0.976, 0.964, 0.973, 0.923},
  {0.945, 0.966, 0.972, 0.977, 0.995, 0.995, 0.986, 1.006, 0.944, 0.977, 1.000, 1.006, 0.996, 0.993, 0.978, 0.966, 0.974, 0.924}};

const double TopElectronSFTool::ele_trigger_errmatrix_e20_medium[6][18] = {
  {0.014, 0.007, 0.007, 0.008, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.014},
  {0.014, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.013},
  {0.014, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.013},
  {0.014, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.013},
  {0.014, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.013},
  {0.014, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.013}};

const double TopElectronSFTool::ele_trigger_errmatrix_e22_medium[6][18] = {
  {1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000},
  {0.018, 0.008, 0.009, 0.009, 0.008, 0.008, 0.008, 0.008, 0.009, 0.014, 0.013, 0.013, 0.013, 0.013, 0.008, 0.009, 0.008, 0.017},
  {0.016, 0.006, 0.006, 0.007, 0.006, 0.006, 0.006, 0.006, 0.007, 0.012, 0.012, 0.012, 0.012, 0.012, 0.006, 0.006, 0.006, 0.016},
  {0.016, 0.006, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.012, 0.012, 0.012, 0.012, 0.012, 0.006, 0.006, 0.006, 0.016},
  {0.016, 0.006, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.012, 0.011, 0.011, 0.011, 0.012, 0.005, 0.006, 0.005, 0.016},
  {0.016, 0.005, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.012, 0.011, 0.011, 0.011, 0.012, 0.005, 0.006, 0.005, 0.016}};

const double TopElectronSFTool::ele_trigger_errmatrix_e22vh_medium[6][18] = {
  {1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000},
  {0.014, 0.007, 0.007, 0.008, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.014},
  {0.013, 0.006, 0.006, 0.007, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.014},
  {0.013, 0.005, 0.005, 0.007, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.005, 0.005, 0.014},
  {0.013, 0.005, 0.005, 0.007, 0.005, 0.005, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.005, 0.005, 0.014},
  {0.013, 0.005, 0.005, 0.007, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.014}};


const double TopElectronSFTool::ele_ID_AFII_SFmatrix[8][18] = {
  {0.949, 0.915, 0.890, 0.902, 0.904, 0.921, 0.923, 0.926, 0.916, 0.944, 0.925, 0.923, 0.923, 0.938, 0.913, 0.921, 0.919, 0.952},
  {0.964, 0.959, 0.947, 0.940, 0.942, 0.955, 0.947, 0.951, 0.993, 0.975, 0.973, 0.957, 0.951, 0.953, 0.935, 0.954, 0.958, 0.971},
  {1.011, 0.983, 0.950, 0.946, 0.974, 0.970, 0.974, 0.981, 0.994, 0.990, 0.979, 0.971, 0.963, 0.964, 0.956, 0.958, 0.981, 1.019},
  {1.029, 1.000, 0.979, 0.951, 0.965, 0.972, 0.977, 0.987, 1.007, 1.018, 0.988, 0.978, 0.966, 0.962, 0.967, 0.960, 0.992, 1.023},
  {1.039, 1.006, 0.964, 0.957, 0.956, 0.974, 0.971, 0.986, 1.017, 1.005, 0.988, 0.974, 0.979, 0.959, 0.966, 0.981, 1.003, 1.032},
  {1.041, 1.014, 0.980, 0.965, 0.962, 0.981, 0.975, 0.990, 1.004, 0.996, 0.992, 0.981, 0.983, 0.963, 0.976, 0.980, 0.997, 1.060},
  {1.052, 1.017, 0.981, 0.959, 0.957, 0.971, 0.973, 0.979, 1.010, 0.997, 0.988, 0.972, 0.969, 0.961, 0.969, 0.979, 0.996, 1.027},
  {1.047, 1.020, 0.985, 0.963, 0.939, 0.969, 0.960, 0.958, 1.032, 1.000, 0.971, 0.961, 0.961, 0.948, 0.965, 1.008, 1.001, 1.028}};

const double TopElectronSFTool::ele_ID_AFII_errmatrix[8][18] = {
  {0.043, 0.040, 0.037, 0.037, 0.036, 0.036, 0.036, 0.035, 0.041, 0.042, 0.036, 0.037, 0.036, 0.036, 0.037, 0.037, 0.040, 0.044},
  {0.037, 0.035, 0.031, 0.031, 0.029, 0.029, 0.029, 0.029, 0.036, 0.036, 0.030, 0.030, 0.029, 0.029, 0.031, 0.031, 0.034, 0.036},
  {0.037, 0.035, 0.031, 0.031, 0.029, 0.029, 0.029, 0.029, 0.036, 0.036, 0.030, 0.030, 0.029, 0.029, 0.031, 0.031, 0.035, 0.037},
  {0.038, 0.035, 0.032, 0.031, 0.030, 0.029, 0.029, 0.029, 0.036, 0.036, 0.030, 0.030, 0.029, 0.030, 0.031, 0.031, 0.035, 0.037},
  {0.037, 0.035, 0.031, 0.031, 0.029, 0.029, 0.029, 0.029, 0.036, 0.036, 0.030, 0.030, 0.029, 0.029, 0.031, 0.031, 0.035, 0.037},
  {0.038, 0.035, 0.032, 0.031, 0.029, 0.029, 0.029, 0.029, 0.036, 0.036, 0.030, 0.030, 0.029, 0.029, 0.031, 0.031, 0.035, 0.039},
  {0.042, 0.036, 0.032, 0.031, 0.030, 0.029, 0.029, 0.029, 0.036, 0.036, 0.030, 0.030, 0.029, 0.030, 0.031, 0.032, 0.036, 0.041},
  {0.053, 0.039, 0.035, 0.033, 0.031, 0.030, 0.031, 0.030, 0.039, 0.039, 0.031, 0.032, 0.030, 0.031, 0.033, 0.035, 0.038, 0.050}};

const double TopElectronSFTool::ele_trigger_AFII_SFmatrix_e20_medium[6][18] = {
  {0.974, 0.975, 0.983, 0.982, 0.996, 0.989, 0.975, 1.002, 0.917, 0.970, 1.004, 1.001, 0.997, 1.004, 0.984, 0.979, 0.983, 0.976},
  {0.967, 0.967, 0.976, 0.975, 0.989, 0.982, 0.968, 0.995, 0.910, 0.963, 0.997, 0.994, 0.989, 0.997, 0.977, 0.972, 0.976, 0.969},
  {0.970, 0.971, 0.979, 0.978, 0.992, 0.985, 0.971, 0.998, 0.913, 0.966, 1.000, 0.997, 0.992, 1.000, 0.980, 0.975, 0.979, 0.972},
  {0.971, 0.972, 0.981, 0.980, 0.994, 0.986, 0.973, 1.000, 0.915, 0.967, 1.002, 0.999, 0.994, 1.001, 0.982, 0.977, 0.980, 0.974},
  {0.974, 0.975, 0.984, 0.983, 0.997, 0.989, 0.976, 1.003, 0.918, 0.970, 1.005, 1.002, 0.997, 1.004, 0.985, 0.980, 0.983, 0.977},
  {0.976, 0.977, 0.986, 0.985, 0.999, 0.992, 0.978, 1.005, 0.920, 0.973, 1.007, 1.004, 0.999, 1.007, 0.987, 0.982, 0.986, 0.979}};

const double TopElectronSFTool::ele_trigger_AFII_SFmatrix_e22_medium[6][18] = {
  {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
  {1.003, 0.974, 0.982, 0.978, 0.993, 0.986, 0.973, 0.998, 0.912, 0.983, 1.018, 1.017, 1.012, 1.014, 0.982, 0.970, 0.982, 0.990},
  {1.001, 0.971, 0.979, 0.975, 0.990, 0.984, 0.971, 0.996, 0.910, 0.981, 1.015, 1.015, 1.009, 1.012, 0.980, 0.968, 0.980, 0.988},
  {1.006, 0.976, 0.984, 0.980, 0.995, 0.988, 0.976, 1.001, 0.915, 0.985, 1.020, 1.020, 1.014, 1.017, 0.984, 0.973, 0.985, 0.993},
  {1.007, 0.977, 0.985, 0.982, 0.996, 0.990, 0.977, 1.002, 0.916, 0.987, 1.022, 1.021, 1.015, 1.018, 0.986, 0.974, 0.986, 0.994},
  {1.009, 0.979, 0.987, 0.983, 0.998, 0.991, 0.979, 1.004, 0.917, 0.988, 1.023, 1.023, 1.017, 1.020, 0.987, 0.976, 0.988, 0.995}};

const double TopElectronSFTool::ele_trigger_AFII_SFmatrix_e22vh_medium[6][18] = {
  {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000},
  {0.924, 0.962, 0.974, 0.979, 0.994, 0.988, 0.974, 1.001, 0.935, 0.966, 0.994, 0.994, 0.988, 0.993, 0.979, 0.967, 0.970, 0.912},
  {0.927, 0.964, 0.976, 0.982, 0.997, 0.991, 0.977, 1.003, 0.938, 0.969, 0.997, 0.997, 0.991, 0.996, 0.981, 0.969, 0.973, 0.914},
  {0.925, 0.963, 0.975, 0.980, 0.995, 0.990, 0.975, 1.002, 0.936, 0.968, 0.995, 0.995, 0.990, 0.994, 0.980, 0.968, 0.971, 0.913},
  {0.927, 0.965, 0.977, 0.982, 0.997, 0.991, 0.977, 1.004, 0.938, 0.969, 0.997, 0.997, 0.991, 0.996, 0.981, 0.969, 0.973, 0.914},
  {0.931, 0.969, 0.981, 0.987, 1.002, 0.996, 0.981, 1.008, 0.942, 0.974, 1.002, 1.001, 0.996, 1.000, 0.986, 0.974, 0.977, 0.919}};

const double TopElectronSFTool::ele_trigger_AFII_errmatrix_e20_medium[6][18] = {
  {0.014, 0.007, 0.007, 0.008, 0.007, 0.007, 0.007, 0.007, 0.007, 0.008, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.013},
  {0.013, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.013},
  {0.013, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.012},
  {0.013, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.012},
  {0.013, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.012},
  {0.013, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.012}};

const double TopElectronSFTool::ele_trigger_AFII_errmatrix_e22_medium[6][18] = {
  {1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000},
  {0.015, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.015},
  {0.015, 0.006, 0.006, 0.007, 0.006, 0.006, 0.006, 0.006, 0.007, 0.007, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.014},
  {0.015, 0.006, 0.006, 0.007, 0.006, 0.006, 0.006, 0.006, 0.007, 0.007, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.014},
  {0.014, 0.005, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.005, 0.014},
  {0.014, 0.005, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.005, 0.014}};

const double TopElectronSFTool::ele_trigger_AFII_errmatrix_e22vh_medium[6][18] = {
  {1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000},
  {0.013, 0.007, 0.007, 0.008, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.014},
  {0.013, 0.006, 0.006, 0.007, 0.006, 0.006, 0.006, 0.006, 0.006, 0.007, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.013},
  {0.013, 0.006, 0.006, 0.006, 0.006, 0.006, 0.005, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.013},
  {0.013, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.013},
  {0.013, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.013}};


// For ID+iso scale factors
double TopElectronSFTool::ele_ID_SF(double eta, double ET)
{

    double ele_eta = eta;
    double ele_ET = ET;
    
    // do range checks first and get out before allocating the arrays ...
    if ( std::fabs(ele_eta) > ele_ID_etamax || (std::fabs(ele_eta) > ele_ID_etacrack[0] && std::fabs(ele_eta) < ele_ID_etacrack[1]) || ele_ET < ele_ID_EtMin) return 1.0;

    int etaI=-1; int ET_I=-1;   
    for (int i=8; i>=0; i--){    // find eta index
        if ( ele_eta > ele_ID_etabins[i] ) {
            etaI = i;
            break;
        }
    }
    for (int i=9; i>=0; i--){    // find eta index
        if ( ele_ET > ele_ID_ETbins[i] ) {
            ET_I = i;
            break;
        }
    }

    if (etaI == -1 || ET_I == -1) return 1.0; // given the range check above this should not happen ...
    return ele_ID_SFmatrix[ET_I][etaI];
}


// For ID+iso scale factor uncertainties (symmetric)
double TopElectronSFTool::ele_ID_SF_err(double eta, double ET)
{

    double ele_eta = eta;
    double ele_ET = ET;
    
    // do range checks first and get out before allocating the arrays ...
    if ( std::fabs(ele_eta) > ele_ID_etamax || (std::fabs(ele_eta) > ele_ID_etacrack[0] && std::fabs(ele_eta) < ele_ID_etacrack[1]) || ele_ET < ele_ID_EtMin) return 1.0;

    int etaI=-1; int ET_I=-1;
    for (int i=8; i>=0; i--){    // find eta index
        if ( ele_eta > ele_ID_etabins[i] ) {
            etaI = i;
            break;
        }
    }
    for (int i=9; i>=0; i--){    // find ET index
        if ( ele_ET > ele_ID_ETbins[i] ) {
            ET_I = i;
            break;
        }
    }

    if (etaI == -1 || ET_I == -1) return 1.0; // given the range check above this should not happen ...
    return ele_ID_errmatrix[ET_I][etaI];
}


double TopElectronSFTool::ele_reco_SF(double eta)
{
    double ele_eta = eta;

    // do range checks first and get out before allocating the arrays ...
    if ( std::fabs(ele_eta) > ele_reco_etamax || (std::fabs(ele_eta) > ele_reco_etacrack[0] && std::fabs(ele_eta) < ele_reco_etacrack[1]) )  return 1.0;

    int etaI=-1;
    for (int i=8; i>=0; i--){    // find eta index
        if ( ele_eta > ele_reco_etabins[i] ) {
            etaI = i;
            break;
        }
    }

    if (etaI == -1) return 1.0; // given the range check above this should not happen ...
    return ele_reco_SFmatrix[etaI];
}


// For reco SF uncertainties (symmetric)
double TopElectronSFTool::ele_reco_SF_err(double eta)
{
    double ele_eta = eta;

    // do range checks first and get out before allocating the arrays ...
    if ( std::fabs(ele_eta) > ele_reco_etamax || (std::fabs(ele_eta) > ele_reco_etacrack[0] && std::fabs(ele_eta) < ele_reco_etacrack[1]) ) return 1.0;

    int etaI=-1; 
    for (int i=8; i>=0; i--){    // find eta index
        if ( ele_eta > ele_reco_etabins[i] ) {
            etaI = i;
            break;
        }
    }

    if (etaI == -1) return 1.0; // given the range check above this should not happen ...
    return ele_reco_errmatrix[etaI];
}


//returns cumulative reco+ID SF
double TopElectronSFTool::ele_recoID_SF(double eta, double ET)
{

	double ele_eta = eta;
	double ele_ET = ET;
	double total_SF = ele_ID_SF(ele_eta, ele_ET) * ele_reco_SF(ele_eta);

	return total_SF;

}

//returns cumulative reco+ID SF uncertainty
double TopElectronSFTool::ele_recoID_SF_err(double eta, double ET)
{

	double ele_eta = eta;
	double ele_ET = ET;
	double ID_err = ele_ID_SF_err(ele_eta, ele_ET) / ele_ID_SF(ele_eta, ele_ET);  //need relative errors
	double reco_err = ele_reco_SF_err(ele_eta) / ele_reco_SF(ele_eta);

	double tot_rel_err = sqrt( pow(ID_err,2) + pow(reco_err,2) );
	double tot_abs_err = tot_rel_err * ele_recoID_SF(ele_eta, ele_ET);

	return tot_abs_err;

}

// For trigger SFs
double TopElectronSFTool::ele_trigger_SF(double eta, double ET, int set)
{
    if(set<0 || set>2) {
      std::cout << "electron_SF_R17.h ele_trigger_SF(): the trigger needs to be 0, 1 or 2. You gave " << set << ". Returning 0. as trigger scale factor." << std::endl;
      return 0.;
    }

    double ele_eta = eta;
    double ele_ET = ET;

    // do range checks first and get out before allocating the arrays ...
    if ( std::fabs(ele_eta) > ele_trigger_etamax || (std::fabs(ele_eta) > ele_trigger_etacrack[0] && std::fabs(ele_eta) < ele_trigger_etacrack[1]) || ele_ET < ele_trigger_EtMin) return 1.0;

    int etaI=-1; int ET_I=-1;
    for (int i=17; i>=0; i--){    // find eta index
        if ( ele_eta > ele_trigger_etabins[i] ) {
            etaI = i;
            break;
        }
    }
    for (int i=5; i>=0; i--){    // find ET index
        if ( ele_ET > ele_trigger_ETbins[i] ) {
            ET_I = i;
            break;
        }
    }
    if (etaI == -1 || ET_I == -1) return 1.0; // given the range check above this should not happen ...

    if(set==0)
      return ele_trigger_SFmatrix_e20_medium[ET_I][etaI];
    else if(set==1)
      return ele_trigger_SFmatrix_e22_medium[ET_I][etaI];
    else if(set==2)
      return ele_trigger_SFmatrix_e22vh_medium[ET_I][etaI];
    else
      return 0.;
}

// For trigger SF uncertainties (symmetric)
double TopElectronSFTool::ele_trigger_SF_err(double eta, double ET, int set)
{
    if(set<0 || set>2) {
      std::cout << "electron_SF_R17.h ele_trigger_SF_err(): the trigger needs to be 0, 1 or 2. You gave " << set << ". Returning 0. as trigger error scale factor." << std::endl;
      return 0.;
    }

    double ele_eta = eta;
    double ele_ET = ET;

    // do range checks first and get out before allocating the arrays ...
    if ( std::fabs(ele_eta) > ele_trigger_etamax || (std::fabs(ele_eta) > ele_trigger_etacrack[0] && std::fabs(ele_eta) < ele_trigger_etacrack[1]) || ele_ET < ele_trigger_EtMin) return 1.0;
    
    int etaI=-1; int ET_I=-1;
    for (int i=17; i>=0; i--){    // find eta index
        if ( ele_eta > ele_trigger_etabins[i] ) {
            etaI = i;
            break;
        }
    }
    for (int i=5; i>=0; i--){    // find ET index
        if ( ele_ET > ele_trigger_ETbins[i] ) {
            ET_I = i;
            break;
        }
    }
    if (etaI == -1 || ET_I == -1) return 1.0; // given the range check above this should not happen ...

    if(set==0)
      return ele_trigger_errmatrix_e20_medium[ET_I][etaI];
    else if(set==1)
      return ele_trigger_errmatrix_e22_medium[ET_I][etaI];
    else if(set==2)
      return ele_trigger_errmatrix_e22vh_medium[ET_I][etaI];
    else
      return 0;
}


// For ID+iso scale factors
double TopElectronSFTool::ele_ID_SF_AFII(double eta, double ET)
{

    double ele_eta = eta;
    double ele_ET = ET;
    int etaI=-1; int ET_I=-1;   

    if ( std::fabs(ele_eta) > ele_ID_etamax || (std::fabs(ele_eta) > ele_ID_etacrack[0] && std::fabs(ele_eta) < ele_ID_etacrack[1]) )  // check forward, crack regions

        return 1.0;       

    else {           
        
        for (int i=17; i>=0; i--){    // find eta index
            if ( ele_eta > ele_ID_etabins[i] ) {
                etaI = i;
                break;
            }
        }
        for (int i=7; i>=0; i--){    // find eta index
            if ( ele_ET > ele_ID_ETbins[i] ) {
                ET_I = i;
                break;
            }
        }

        return ele_ID_AFII_SFmatrix[ET_I][etaI];

    } //else

    return 0;
}

// For ID+iso scale factor uncertainties (symmetric)
double TopElectronSFTool::ele_ID_SF_err_AFII(double eta, double ET)
{

    double ele_eta = eta;
    double ele_ET = ET;
    int etaI=-1; int ET_I=-1;   

    if ( std::fabs(ele_eta) > ele_ID_etamax || (std::fabs(ele_eta) > ele_ID_etacrack[0] && std::fabs(ele_eta) < ele_ID_etacrack[1]) )  // check forward, crack regions

        return 1.0;       

    else {           
        
        for (int i=17; i>=0; i--){    // find eta index
            if ( ele_eta > ele_ID_etabins[i] ) {
                etaI = i;
                break;
            }
        }
        for (int i=7; i>=0; i--){    // find ET index
            if ( ele_ET > ele_ID_ETbins[i] ) {
                ET_I = i;
                break;
            }
        }
		
        return ele_ID_AFII_errmatrix[ET_I][etaI];

    } //else

    return 0;
}


//returns cumulative reco+ID SF
double TopElectronSFTool::ele_recoID_SF_AFII(double eta, double ET)
{

	double ele_eta = eta;
	double ele_ET = ET;
	double total_SF = ele_ID_SF_AFII(ele_eta, ele_ET) * ele_reco_SF(ele_eta);

	return total_SF;

}

//returns cumulative reco+ID SF uncertainty
double TopElectronSFTool::ele_recoID_SF_err_AFII(double eta, double ET)
{

	double ele_eta = eta;
	double ele_ET = ET;
	double ID_err = ele_ID_SF_err_AFII(ele_eta, ele_ET) / ele_ID_SF_AFII(ele_eta, ele_ET);  //need relative errors
	double reco_err = ele_reco_SF_err(ele_eta) / ele_reco_SF(ele_eta);

	double tot_rel_err = sqrt( pow(ID_err,2) + pow(reco_err,2) );
	double tot_abs_err = tot_rel_err * ele_recoID_SF(ele_eta, ele_ET);

	return tot_abs_err;

}

// For trigger SFs
double TopElectronSFTool::ele_trigger_SF_AFII(double eta, double ET, int set)
{

    double ele_eta = eta;
    double ele_ET = ET;
    int etaI=-1; int ET_I=-1;   

    if ( std::fabs(ele_eta) > ele_trigger_etamax || (std::fabs(ele_eta) > ele_trigger_etacrack[0] && std::fabs(ele_eta) < ele_trigger_etacrack[1]) )  // check forward, crack regions

        return 1.0;       

    else {           
        
        for (int i=17; i>=0; i--){    // find eta index
            if ( ele_eta > ele_trigger_etabins[i] ) {
                etaI = i;
                break;
            }
        }
        for (int i=5; i>=0; i--){    // find ET index
            if ( ele_ET > ele_trigger_ETbins[i] ) {
                ET_I = i;
                break;
            }
        }
		
       	if(set==0)
	  return ele_trigger_AFII_SFmatrix_e20_medium[ET_I][etaI];
	else if(set==1)
	  return ele_trigger_AFII_SFmatrix_e22_medium[ET_I][etaI];
	else if(set==2)
	  return ele_trigger_AFII_SFmatrix_e22vh_medium[ET_I][etaI];
	else 
	  return 0;
	
    } //else

    return 0;

}

// For trigger SF uncertainties (symmetric)
double TopElectronSFTool::ele_trigger_SF_err_AFII(double eta, double ET, int set)
{

    double ele_eta = eta;
    double ele_ET = ET;
    int etaI=-1; int ET_I=-1;   

    if ( std::fabs(ele_eta) > ele_trigger_etamax || (std::fabs(ele_eta) > ele_trigger_etacrack[0] && std::fabs(ele_eta) < ele_trigger_etacrack[1]) )  // check forward, crack regions

        return 1.0;       

    else {           
        
        for (int i=17; i>=0; i--){    // find eta index
            if ( ele_eta > ele_trigger_etabins[i] ) {
                etaI = i;
                break;
            }
        }
        for (int i=5; i>=0; i--){    // find ET index
            if ( ele_ET > ele_trigger_ETbins[i] ) {
                ET_I = i;
                break;
            }
        }
	
	if(set==0)
	  return ele_trigger_AFII_errmatrix_e20_medium[ET_I][etaI];
	else if(set==1)
	  return ele_trigger_AFII_errmatrix_e22_medium[ET_I][etaI];
	else if(set==2)
	  return ele_trigger_AFII_errmatrix_e22vh_medium[ET_I][etaI];
	else 
	  return 0;
    } //else

    return 0;
}
