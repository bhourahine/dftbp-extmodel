/**
   @brief
   Headers for an example external model for the DFTB+ project:
   www.dftbplus.org Copyright (C) 2022 B. Hourahine

   See the LICENSE file for terms of usage and distribution.

   @file
*/

#ifndef H_EXTERNALMODEL
#define H_EXTERNALMODEL

#ifdef __cplusplus
extern "C" {
#endif

  /** Declared capabilities of this model for initial setup of the
      DFTB+ code */
  typedef struct {

    /** Can this library deliver a hamiltonian to the calling code? */
    bool hamiltonian;

    /** Can this library deliver an overlap matrix (i.e. is a
	non-orthogonal model)? */
    bool overlap;

    /** Does this library deliver energy terms (other than the
	band-structure energy)? */
    bool energy;

    /** Order of derivatives returned by the model */
    int derivativeOrder;

    /** Is the model self-consistent */
    bool selfconsistent;

    /** Number of spin channels supported by model (none:0,
	collinear:1, non-collinear:3) */
    int spinchannels;

  } mycapabilities;


/** Internal state of this model, including any initialisation passed
    in from DFTB+, any data read in, checks for state, etc.
*/
struct mystate {

  // Is this structure initialised
  bool initialised;

  // References to data structures from DFTB+

  // Chemical species of atoms in global geometric structure
  int* globalSpeciesOfAtoms;

  // Clusters around atom sites

  // Number of atom centred clusters
  int nAtomicClusters;
  // Index for where the cluster starts in the following two arrays
  int* indexAtomicClusters;
  // Geometries of clusters
  double* atomicClusters;
  // Numbering for the atoms of the cluster in the original (global)
  // geometry of the system
  int* atomicGlobalAtNos;
  // Indexing for diatomic elements of H or S matrices
  int* atomClusterIndex;

  // Clusters around diatomic bonds

  // Number of bond clusters
  int nBndClusters;
  // Index for where the cluster starts in the following two arrays
  int* indexBndClusters;
  // coordinates in the clusters
  double* bndClusters;
  // Numbering for the atoms of the cluster in the original (global)
  // geometry of the system
  int* bndGlobalAtNos;
  // Indexing for diatomic elements of H or S matrices
  int* bondClusterIndex;

  // internal model parameters
  double onsites[3]; // Hs, Cs and Cs* alpha
  double hopping[6]; // Hs-Hs, Hs-Cs, Hs-Cs*, Cs-Cs, Cs-Cs*, Cs*-Cs*
		     // beta values
  // index from order of chemical species in DFTB+ to this code's
  // parameter ordering:
  int species2params[2];

};


  /**
     Combined application programming interface and application binary
     interface semantic version, as supported by this library. Note,
     the format and is externally defined in the DFTB+ project

     @param major Major version, revised on breaking changes
     @param minor Minor version, revised on extensions
     @param patch Patch version, revised on invisible changes

  */
  void dftbp_model_apbi(int* major, int* minor, int* patch);

  /**
     Declare capabilities of this model to DFTB+ via the external model API.

     @param modelname null terminated string for name of this model
     @param nChars characters in model name
     @param capabilities structure with capabilities of the model

  */
  int dftbp_provided_with(typeof (mycapabilities) *capabilities, int* nChars, char** modelname);

  /**
     Set up this model, read some settings from DFTB+ over it's
     external model API and initialise it's data structure for
     handling via DFTB+.

     @param nspecies number of chemical species/types present

     @param species array of null terminated strings labelling
     chemical species

     @param interactionCutoff Longest cutoff for, i.e. distance over
     which atoms of each species have hamiltonian or repulsive
     interactions

     @param environmentCutoff Distance over which neighbours influence
     interactions, i,e. 0 if model is environmentally independent, or
     nearest-neighbour/longer if surrounding atoms influence
     interactions. This is used to cut out local clusters of around
     the interacting atom/dimer

     @param nShellsOnSpecies number of shells of atomic orbitals, set
     to 0 if not a hamiltonian model

     @param shells Angular momentum of shells species resolved atomic
     shells, freed on return to DFTB+

     @param shellOccs Reference occupation for neutrality, freed on
     return to DFTB+

     @param state internal state and data of the model, not checked by
     DFTB+, just passed around

     @param message return message, in event of routine failure
     (return != 0)

     @return 0 on successful return, non-zero if there is an error
     message to check

   */
  int initialise_model_for_dftbp(int* nspecies, char* speciesName[],
				 double* interactionCutoff,
				 double* environmentCutoff,
				 int* nShellsOnSpecies[], int** shells,
                                 double** shellOccs, intptr_t *state,
				 char** message);


    /**
     Update this model, using geometric and other information from
     DFTB+ over it's external model API.

     @param state internal state and data of the model, this is not
     checke by DFTB+, just passed around by it

     @param species Species index for atoms in the global structure

     @param nAtomicClusters Number of atom centred clusters

     @param indexAtomicClusters starting index for location of
     coordinates in the atomicClusters array

     @param atomicClusters Geometric clusters centred on atoms for the
     onsite matrix element predictions

     @param atomicGlobalAtNos Numbers of atoms from clusters in the
     global system

     @param nBndClusters Number of bond centred clusters

     @param indexBndClusters starting index for location of
     coordinates in the bndClusters array

     @param bndClusters Geometric clusters centred on bonds for the
     diatomic matrix element predictions

     @param bndGlobalAtNos Numbers of atoms from bond clusters in the
     global system

     @param atomClusterIndex Indexing for where to find atomic
     (onsite) elements in hamiltonian/overlap matrices

     @param bondClusterIndex Indexing for where to find diatomic
     (bond) elements in hamiltonian/overlap matrices

     @param message return message, in event of routine failure
     (return != 0)

     @return 0 on successful return, non-zero if there is an error
     message to check

   */
  int update_model_for_dftbp(intptr_t *state, int* species,
                             int* nAtomicClusters, int* indexAtomicClusters,
                             double* atomicClusters, int* atomicGlobalAtNos,
                             int* nBndClusters, int* indexBndClusters,
			     double* bndClusters, int* bndGlobalAtNos,
			     int* atomClusterIndex, int* bondClusterIndex,
                             char** message);


  /**
      Get model predictions

      @param state internal state and data of the model, this is not
      checke by DFTB+, just passed around by it

      @param h0 hamiltonian

      @param over overlap matrix

      @param h0Index hamiltonian index for blocks of matrix elements

      @param message return message, in event of routine failure
      (return != 0)

      @return 0 on successful return, non-zero if there is an error
      message to check

  */
  int predict_model_for_dftbp(intptr_t *state, double *h0, double *over,
                              char** message);


  /**
     Clean up after this model, freeing any memory in the mystate type

     @param state internal state and data of the model. This is not
     checke by DFTB+, just passed around by it, so we need to remove
     any allocated memory here.

     @param message return message, in event of routine failure
     (return != 0)

     @return 0 on successful return, non-zero if there is an error
     message to check

  */
  void cleanup_model_for_dftbp(intptr_t *state);
  
#ifdef __cplusplus
}
#endif

#endif
