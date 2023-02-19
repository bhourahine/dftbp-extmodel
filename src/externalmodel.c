/**

@brief
Example external model for the DFTB+ project: www.dftbplus.org
Copyright (C) 2022 - 2023 B. Hourahine

See the LICENSE file for terms of usage and distribution.

@file
*/

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/asprintf.h" // As this function isn't cross platform available
#include "../include/externalmodel.h" // Structures for DFTB+ communication

// Project version
void dftbp_model_apbi(int* major, int* minor, int* patch)
{
  *major = 0;
  *minor = 1;
  *patch = 0;
}


// Declare capabilities of this model to DFTB+
int dftbp_provided_with(typeof (mycapabilities) *capabilities, int* nChars, char** modelname)
{

  // flags for capabilities of this specific model
  *capabilities = (mycapabilities) {
    .hamiltonian = true, .overlap = false, .energy = false,
    .derivativeOrder = 0, .selfconsistent = false, .spinchannels = 0
  };

  int lenName = asprintf(modelname, "%s", "Huckel toy model");

  if (lenName != -1) {
    *nChars = lenName;
    return 1;
  } else {
    *nChars = 0;
    // bring down the code messily
    return -1;
  }

}


// Initialise this model, reading some settings from DFTB+
int initialise_model_for_dftbp(int* nspecies, char* speciesName[],
                               double* interactionCutoff, double* environmentCutoff,
                               int** nShellsOnSpecies, int** shellLValues, double** shellOccs,
                               intptr_t *state, int* nChars, char** message)
{

  // Allocate structure for internal state, and generate am intptr for
  // return to DFTB+
  struct mystate* internalState = (struct mystate*) malloc(sizeof(struct mystate));
  *state = (intptr_t) internalState;

  FILE *input;
  int ii, items, iStat;

  // Open input file for some constants for this model, assuming it's
  // in the runtime directory
  input = fopen("input.dat", "r");
  if (!input) {
    iStat = asprintf(message, "%s", "Library error opening input file");
    if (iStat != -1) {
      *nChars = iStat;
    } else {
      *nChars = 0;
    }
    // bring down the code messily
    return -1;
  }

  // read ancillary input file for model parameters and then store
  // them into the model's internal structure, that in turn will get passed
  // around between calls to the model.
  items = fscanf(input, "%lf %lf %lf", &internalState->onsites[0],
                 &internalState->onsites[1], &internalState->onsites[2]);
  if (items == EOF) {
    iStat = asprintf(message, "%s", "Toy library malformed end of data file at first line");
    if (iStat != -1) {
      *nChars = iStat;
    } else {
      *nChars = 0;
    }
    // bring down the code messily
    return -1;
  }
  if (items != 3) {
    iStat = asprintf(message, "%s %i", "Toy library malformed first line of data file:", items);
    if (iStat != -1) {
      *nChars = iStat;
    } else {
      *nChars = 0;
    }
    // bring down the code messily
    return -1;
  }
  items = fscanf(input, "%lf %lf %lf %lf %lf %lf", &internalState->hopping[0],
                 &internalState->hopping[1],
                 &internalState->hopping[2],
                 &internalState->hopping[3],
                 &internalState->hopping[4],
                 &internalState->hopping[5]);
  if (items == EOF) {
    iStat = asprintf(message,"%s", "Toy library malformed end of data file before 2nd line");
    if (iStat != -1) {
      *nChars = iStat;
    } else {
      *nChars = 0;
    }
    // bring down the code messily
    return -1;
  }
  if (items != 6) {
    iStat = asprintf(message, "%s", "Toy library malformed second line of data file");
    if (iStat != -1) {
      *nChars = iStat;
    } else {
      *nChars = 0;
    }
    // bring down the code messily
    return -1;
  }
  items = fscanf(input, "%lf", interactionCutoff);
  if (items == EOF) {
    iStat = asprintf(message, "%s", "Toy library malformed end of data file before 3rd line (bond cut-off)\n");
    if (iStat != -1) {
      *nChars = iStat;
    } else {
      *nChars = 0;
    }
    // bring down the code messily
    return -1;
  }
  if (items != 1) {
    iStat = asprintf(message, "%s", "Toy library malformed third line of data file\n");
    if (iStat != -1) {
      *nChars = iStat;
    } else {
      *nChars = 0;
    }
    // bring down the code messily
    return -1;
  }

  for (ii = 0; ii < *nspecies; ii++) {
    // This specific model is only for H and C atoms, so will throw an
    // error otherwise
    if (strcmp(speciesName[ii], "C") != 0 && strcmp(speciesName[ii], "H") != 0)
      {
	iStat = asprintf(message, "%s %s", "Toy library only knows about C and H atoms, not",
                        speciesName[ii]);
        if (iStat != -1) {
          *nChars = iStat;
        } else {
          *nChars = 0;
        }
        // bring down the code messily
        return -1;
      }
    if (strcmp(speciesName[ii], "H") == 0) {
      (*internalState).species2params[ii] = 0;
    }
    if (strcmp(speciesName[ii], "C") == 0) {
      (*internalState).species2params[ii] = 1;
    }
  }

  // Surroundings required around atoms and bonds:
  *environmentCutoff = 4.0;

  // Maximum number of shells of orbitals on atoms in this model
  int maxShells = 2;

  /* This particular model is so only a single s shell on H and two s
     shells on C Note count for shells uses Fortran convention
     starting at 1, while L uses physics convention of s=0 */
  *nShellsOnSpecies =  calloc(*nspecies, sizeof(int));
  // Note, this will be a fortran array on the other end, so
  // [nSpecies][maxShells] (hence the 2):
  *shellLValues =  calloc(*nspecies * maxShells, sizeof(int));
  for (ii=0; ii<*nspecies; ii++) {
    if ((*internalState).species2params[ii] == 0)
      {
        // Hydrogen, s
        (*nShellsOnSpecies)[ii] = 1;
        (*shellLValues)[ii] = 0;
      } else
      {
        // Carbon, s and s* shells
        (*nShellsOnSpecies)[ii] = 2;
        (*shellLValues)[ii] = 0;
        (*shellLValues)[ii+1] = 0;
      }
  }

  // each atom is neutral when it's single shell containing only one
  // electron:
  *shellOccs =  calloc(*nspecies * maxShells, sizeof(double));
  for (ii=0; ii<*nspecies; ii++) {
    int jj = ii * maxShells;
    if ((*internalState).species2params[ii] == 0)
      {
        // H atom, the s orbital is occupied
        (*shellOccs)[jj] = 1.0;
      } else {
      // C atom, only the s orbital is occupied
      (*shellOccs)[jj] = 1.0;
    }
  }

  (*internalState).initialised = true;

  (*internalState).nAtomicClusters = 0;
  (*internalState).indexAtomicClusters = NULL;
  (*internalState).atomicClusters = NULL;
  (*internalState).atomicGlobalAtNos = NULL;

  (*internalState).nBndClusters = 0;
  (*internalState).indexBndClusters = NULL;
  (*internalState).bndClusters = NULL;
  (*internalState).bndGlobalAtNos = NULL;
  (*internalState).bondClusterIndex = NULL;

  iStat = asprintf(message,"Model initialised, on-site energies : H %f, C %f %f\n",
                  (*internalState).onsites[0], (*internalState).onsites[1],
                  (*internalState).onsites[2]);
  if (iStat != -1) {
    *nChars = iStat;
    // return that there is a message
    return 1;
  } else {
    *nChars = 0;
    // bring down the code messily
    return -1;
  }

}


// Update this model, using geometric and other information from DFTB+
int update_model_for_dftbp(intptr_t *state, int* species, int* nAtomicClusters,
                           int* indexAtomicClusters, double* atomicClusters,
                           int* atomicGlobalAtNos, int* nBndClusters,
                           int* indexBndClusters, double* bndClusters,
                           int* bndGlobalAtNos, int* atomClusterIndex,
                           int* bondClusterIndex, int* nChars, char** message)
{

  // map pointer back to structure
  struct mystate* internalState = (struct mystate*) *state;

  int iErr;

  if (!(*internalState).initialised) {
    iErr = asprintf(message, "%s", "External model is not properly initialised");
    if (iErr != -1) {
      *nChars = iErr;
      return -1;
    } else {
      *nChars = 0;
      // bring down the code messily
      return -1;
    }
  }

  internalState->globalSpeciesOfAtoms = species;

  internalState->nAtomicClusters = *nAtomicClusters;
  internalState->indexAtomicClusters = indexAtomicClusters;
  internalState->atomicClusters = atomicClusters;
  internalState->atomicGlobalAtNos = atomicGlobalAtNos;
  internalState->atomClusterIndex = atomClusterIndex;

  //printf("Number of atomic clusters: %i\n", *nAtomicClusters);

  internalState->nBndClusters = *nBndClusters;
  internalState->indexBndClusters = indexBndClusters;
  internalState->bndClusters = bndClusters;
  internalState->bndGlobalAtNos = bndGlobalAtNos;
  internalState->bondClusterIndex = bondClusterIndex;

  //printf("Number of bond clusters: %i\n", *nBndClusters);

  /*printf(" Initial on-site energies : H %f, C %f %f\n",
         (*internalState).onsites[0], (*internalState).onsites[1],
         (*internalState).onsites[2]);
  */

  // blank return message if nothing happening
  *nChars = 0;

  // All OK
  return 0;

}


// Make and then get model predictions back to DFTB+
int predict_model_for_dftbp(intptr_t *state, double *h0, double *over,
                            int* nChars, char** message)
{

  struct mystate* internalState = (struct mystate*) *state;

  int iErr;

  if (!(*internalState).initialised) {
    iErr = asprintf(message, "%s", "External model is not properly initialised");
    if (iErr != -1) {
      *nChars = iErr;
      return -1;
    } else {
      *nChars = 0;
      // bring down the code messily
      return -1;
    }
  }

  // For this specific model, overlap is ignored, but would be indexed
  // in the same way as the hamiltonian

  // global atom number == atomic cluster number
  for (int iAt=0; iAt<(*internalState).nAtomicClusters; iAt++) {
    // First atom location in species etc. arrays
    int iStart = (*internalState).atomClusterIndex[iAt];
    int iSpecies = (*internalState).globalSpeciesOfAtoms[iAt];

    // index for this model's parameter storage, vs the atomic species
    // indices from DFTB+ :
    int iParam = (*internalState).species2params[iSpecies-1];

    // simple use of onsite energies
    if (iParam == 0) // H atom
      {
        h0[iStart] = (*internalState).onsites[0]; // s orbital
      } else // C atom
      {
        // 2x2 symmetric matrix of s and s* onsites, reshaped as a 4
        // element vector:
        h0[iStart] = (*internalState).onsites[1]; // s orbital
        h0[iStart+3] = (*internalState).onsites[2]; // s* orbital
      }


    // Toy crystal field model, demonstrating getting positions and
    // species of atoms in the halo around each atomic site in the
    // geometry
    //
    // Model affects s and s* mixing on C atoms, but only due to the
    // presence of any surrounding H atoms within the environment
    // cutoff distance (ignoring any surrounding C). H atoms
    // themselves would be unaffected by crystal field terms from
    // their surroundings, as there is only a single s orbital on
    // those atoms.

    if (iParam == 1) { // there is a C atom at the cluster centre
      // so check it's neighbours.
      // excludes first atom in cluster:
      int jStart = (*internalState).indexAtomicClusters[iAt];
      int jEnd = (*internalState).indexAtomicClusters[iAt+1] - 1;
      //printf("Atom range in cluster %i : < %i\n", jStart, jEnd);
      for (int iNeighbour = jStart; iNeighbour < jEnd; iNeighbour++) {
        int jAt = (*internalState).atomicGlobalAtNos[iNeighbour];
        int jSpecies = (*internalState).globalSpeciesOfAtoms[jAt-1];
        int jParam = (*internalState).species2params[jSpecies-1];
        if (jParam == 0) { // H atom is in suroundings
          int jCoordStart = 3 * iNeighbour; // stride in coordinates
          // distance squared
          double d2 =
            ((*internalState).atomicClusters[jCoordStart]
             *(*internalState).atomicClusters[jCoordStart]
             +(*internalState).atomicClusters[jCoordStart+1]
             *(*internalState).atomicClusters[jCoordStart+1]
             +(*internalState).atomicClusters[jCoordStart+2]
             *(*internalState).atomicClusters[jCoordStart+2]
             );
          // scale such that final term is ~1E-12 by 4.0 a.u. cutoff
          d2 *= -1.0*(6.0/4.0)*(6.0/4.0);
          h0[iStart+1] += 1000.0*exp(d2);
          h0[iStart+2] += 1000.0*exp(d2);
        }
      }
    }

  }

  // off-site (diatomic) elements
  for (int iClust = 0; iClust < (*internalState).nBndClusters; iClust++) {

    // Index for first atom in cluster
    int iStart = (*internalState).indexBndClusters[iClust] - 1;
    // Index for last atom in  cluster
    int iEnd = (*internalState).indexBndClusters[iClust+1] - 1;

    // First two atoms in the cluster, these having the bond between
    // their orbitals
    int iAt = (*internalState).bndGlobalAtNos[iStart] - 1;
    int jAt = (*internalState).bndGlobalAtNos[iStart+1] - 1;
    int iSp = (*internalState).globalSpeciesOfAtoms[iAt];
    int jSp = (*internalState).globalSpeciesOfAtoms[jAt];

    // vector from first to second atom in the bond
    double v[3];
    for (int jj=0; jj < 3; jj++) {
      v[jj] = (*internalState).bndClusters[3*(iStart+1)+jj];
      v[jj] -= (*internalState).bndClusters[3*iStart+jj];
    }
    double d = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

    int iParam = (*internalState).species2params[iSp-1];
    int jParam = (*internalState).species2params[jSp-1];

    if (iParam == jParam) {
      // Homo-nuclear bond

      if (iParam == 0) {
        // H-H bond
	h0[(*internalState).bondClusterIndex[iClust]] =
	  (*internalState).hopping[0];
      } else {
        // C-C bond
	h0[(*internalState).bondClusterIndex[iClust]] =
	  (*internalState).hopping[3];
	h0[(*internalState).bondClusterIndex[iClust]+1] =
	  (*internalState).hopping[4];
	h0[(*internalState).bondClusterIndex[iClust]+2] =
	  (*internalState).hopping[4];
	h0[(*internalState).bondClusterIndex[iClust]+3] =
	  (*internalState).hopping[5];
      }

    } else {
      // Hetero-nuclear bond

      if (iParam == 0) {
        // H-C bond, as iParam != 0
	h0[(*internalState).bondClusterIndex[iClust]] =
	  (*internalState).hopping[1];
	h0[(*internalState).bondClusterIndex[iClust]+1] =
	  (*internalState).hopping[2];
      } else {
        // C-H bond
	h0[(*internalState).bondClusterIndex[iClust]] =
	  (*internalState).hopping[1];
	h0[(*internalState).bondClusterIndex[iClust]+1] =
	  (*internalState).hopping[2];
      }

    }


    //h0[(*internalState).bondClusterIndex[iClust]] = 0.1 * (double) (iClust+1);

    // Atoms in the rest of the bond cluster, if needed
    /*for (int jj = iStart+2; jj < iEnd; jj++) {
      int kAt = (*internalState).bndGlobalAtNos[jj] - 1;
        printf("Atom surrounding bond number %i: %i\n", iClust, kAt);
    }*/

  }

  // blank return message if nothing failing
  *nChars = 0;

  // Everything OK:
  return 0;

}


// Clean up after this model, freeing any memory in the mystate type
int cleanup_model_for_dftbp(intptr_t *state, int* nChars, char** message) {

  // DFTB+ only sees integer pointer "state", so need original
  // structure to clean up
  struct mystate* internalState = (struct mystate*) *state;

  free(internalState);
  *state = (intptr_t) internalState;

  int iStat = asprintf(message,"%s", "Cleaned up model");
  if (iStat != -1) {
    *nChars = iStat;
    // return that there is a message to check
    return 1;
  } else {
    *nChars = 0;
    // bring down the code messily
    return -1;
  }

}
