/** @file
 * @author Pierre-Antoine Morin
 * @brief Fichier d'entête déclarant la fonction @ref branchAndBound()
 */
#ifndef BRANCH_AND_BOUND_H_INCLUDED
/** Chien de garde pour l'inclusion du fichier @ref branch_and_bound.h */
#define BRANCH_AND_BOUND_H_INCLUDED

//#include "wrapper.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/** Fonction implémentant l'algorithme Branch @& Bound
 * @param [in] dataLocation L'adresse d'un objet de la classe @ref Data, dans lequel lire les durées de fabrication des jobs sur les machines
 * @param [out] sequenceLocation L'adresse d'un objet de la classe @ref Sequence, dans lequel écrire l'ordonnancement optimal trouvé
 * @pre La séquence passée en paramètre doit être vide.
 * @return La valeur minimale du critère
 * @see @ref wrapper.h (Fichier d'entête définissant des fonctions permettant d'utiliser du code C++ dans du code C)
 */
unsigned int branchAndBound(const int m, const int n, const long* p);

#ifdef __cplusplus
}
/* extern "C" */
#endif /* __cplusplus */

#endif /* BRANCH_AND_BOUND_H_INCLUDED */
