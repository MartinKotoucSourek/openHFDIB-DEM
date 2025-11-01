/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "dlvoClass.H"

#include "cyclicPlaneInfo.H"

#include <vector>

using namespace Foam;

//---------------------------------------------------------------------------//
dlvoClass::dlvoClass(const scalar& cutOff):
    cutOff_(cutOff, cutOff, cutOff)
{}

dlvoClass::~dlvoClass()
{}

void dlvoClass::setBBoxes(List<std::shared_ptr<boundBox>> bBox)
{
    HashTable<List<vector>,string,Hash<string>> const& cyclicPlaneInfoList = cyclicPlaneInfo::getCyclicPlaneInfo();
    std::vector<vector> cyclicVectors;
    for (auto patch : cyclicPlaneInfoList.toc())
    {
        auto cyclicVector = cyclicPlaneInfo::getCyclicTransVec(patch);
        cyclicVectors.push_back(cyclicVector);
    }

    auto tmpBBoxes = bBox;
    for (auto bboxI : tmpBBoxes)
    {
        for (auto cyclicVector : cyclicVectors)
        {
            auto newBoundBox = std::make_shared<boundBox>();
            newBoundBox->min() = bboxI->min() - cyclicVector;
            newBoundBox->max() = bboxI->max() - cyclicVector;

            auto iter = std::find_if(bBox.begin(), bBox.end(), [this, newBoundBox](std::shared_ptr<boundBox> cBBox) {
                if ((cBBox->min() > newBoundBox->min() && cBBox->min() < newBoundBox->max()) ||
                    (cBBox->max() > newBoundBox->min() && cBBox->max() < newBoundBox->max()))
                {
                    return true;
                }

                if (mag(cBBox->min() - newBoundBox->min()) < cutOff_[0] ||
                    mag(cBBox->min() - newBoundBox->max()) < cutOff_[0])
                {
                    return true;
                }

                if (mag(cBBox->max() - newBoundBox->min()) < cutOff_[0] ||
                    mag(cBBox->max() - newBoundBox->max()) < cutOff_[0])
                {
                    return true;
                }

                return false;
            });

            if (iter == bBox.end())
            {
                tmpBBoxes.append(newBoundBox);
            }
        }
    }

    if (tmpBBoxes.size() != bBox_.size())
    {
        bBox_.clear();
        for (int i = 0; i < tmpBBoxes.size(); i++)
        {
            bBox_.append(std::make_shared<boundBox>());
        }
    }

    for (int i = 0; i < tmpBBoxes.size(); i++)
    {
        bBox_[i]->min() = tmpBBoxes[i]->min() - cutOff_ / 2;
        bBox_[i]->max() = tmpBBoxes[i]->max() + cutOff_ / 2;
    }
}

List<std::shared_ptr<boundBox>> dlvoClass::getBBoxes()
{
    return bBox_;
}

//---------------------------------------------------------------------------//
