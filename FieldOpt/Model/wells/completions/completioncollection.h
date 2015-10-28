/******************************************************************************
 *
 * completioncollection.h
 *
 * Created: 24.09.2015 2015 by einar
 *
 * This file is part of the FieldOpt project.
 *
 * Copyright (C) 2015-2015 Einar J.M. Baumann <einar.baumann@ntnu.no>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 *****************************************************************************/

#ifndef COMPLETIONCOLLECTION_H
#define COMPLETIONCOLLECTION_H

#include <QList>

#include "perforation.h"
#include "inflowcontroldevice.h"

namespace Model {
namespace Wells {
namespace Completions {

/*!
 * \brief The CompletionCollection class collects all completions for the model and checks whether
 * new completions are valid when they are added (e.g. that they don't overlap).
 */
class CompletionCollection
{
public:
    CompletionCollection();

    void AddPerforation(Perforation *perf); //!< Add a new perforation.
    void AddInflowControlDevice(InflowControlDevice *icd); //!< Add a new ICD.

private:
    QList<Perforation *> *perforations_;
    QList<InflowControlDevice *> *inflow_control_devices_;

    bool NewCompletionIsValid(Perforation *perf);
    bool NewCompletionIsValid(InflowControlDevice *icd);
};

}
}
}

#endif // COMPLETIONCOLLECTION_H
