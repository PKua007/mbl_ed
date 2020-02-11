//
// Created by Piotr Kubala on 10/02/2020.
//

#ifndef MBL_ED_CAVITYCONSTANTSREADER_H
#define MBL_ED_CAVITYCONSTANTSREADER_H

#include <istream>

#include "simulation/CavityConstants.h"

/**
 * @brief A class, which reads CavityConstants from a dat file in the format shown below.
 * @details
 * <pre>
 * phi0_value1 cos_for_site_1 wannier_for_site_1 y_for_site_1
 * phi0_value1 cos_for_site_2 wannier_for_site_2 y_for_site_2
 * phi0_value1 cos_for_site_3 wannier_for_site_3 y_for_site_3
 * phi0_value2 cos_for_site_1 wannier_for_site_1 y_for_site_1
 * phi0_value2 cos_for_site_2 wannier_for_site_2 y_for_site_2
 * phi0_value2 cos_for_site_3 wannier_for_site_3 y_for_site_3
 * </pre>
 * Note that multiple repeated @a phi0 values correspond to data for subsequent sites. The number of sites is
 * determined by counting how much first rows have the same phi0. The pattern must be repeated for all (@a phi0)
 * realisations.
 */
class CavityConstantsReader {
private:
    struct Row;

    static std::vector<Row> loadRows(std::istream &in);
    static std::size_t countNumberOfSites(const std::vector<Row> &rows);

public:
    static CavityConstants load(std::istream &in);
};


#endif //MBL_ED_CAVITYCONSTANTSREADER_H
