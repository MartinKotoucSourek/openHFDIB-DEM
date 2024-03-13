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
#include "vCntcList.H"

using namespace Foam;

//---------------------------------------------------------------------------//
vCntcList::vCntcList()
:
verletList()
{}
vCntcList::~vCntcList()
{}
//---------------------------------------------------------------------------//
void vCntcList::addBodyToVList(immersedBody& ib)
{
    List<std::shared_ptr<boundBox>> bBoxes = ib.getGeomModel().getBBoxes();

    forAll(bBoxes, bBox)
    {
        verletBoxes_.push_back(verletBox::create(
            ib.getBodyId(),
            bBoxes[bBox],
            ib.getbodyOperation()==0
        ));

        verletBoxes_.back()->setVerletPoints();

        forAll(verletLists_, vListI)
        {
            if(vListI != emptyDim)
            {
                verletLists_[vListI].push_back(
                    verletBoxes_.back()->getMinPoint()
                );

                verletLists_[vListI].push_back(
                    verletBoxes_.back()->getMaxPoint()
                );
            }
        }
    }
}
//---------------------------------------------------------------------------//
// void vCntcList::removeBodyFromVList(immersedBody& ib)
// {
//     forAll (verletLists_, coordI)
//     {
//         verletLists_[coordI].remove_if(
//             [&ib](std::shared_ptr<verletPoint>& vPoint)
//             {
//                 return vPoint->getBodyId() == ib.getBodyId();
//             }
//         );

//         for (auto it = cntNeighList_[coordI].begin();
//             it != cntNeighList_[coordI].end();)
//         {
//             if (it->first.first == ib.getBodyId()
//                 ||
//                 it->first.second == ib.getBodyId())
//             {
//                 it = cntNeighList_[coordI].erase(it);
//             }
//             else
//             {
//                 ++it;
//             }
//         }
//     }

//     verletBoxes_.remove_if(
//         [&ib](std::shared_ptr<verletBox>& vBox)
//         {
//             return vBox->getBodyId() == ib.getBodyId();
//         }
//     );

//     for (auto iter = posCntList_.begin(); iter != posCntList_.end();)
//     {
//         if (iter->first == ib.getBodyId()
//             ||
//             iter->second == ib.getBodyId())
//         {
//             iter = posCntList_.erase(iter);
//         }
//         else
//         {
//             ++iter;
//         }
//     }
// }
//---------------------------------------------------------------------------//
void vCntcList::update(PtrList<immersedBody>& ibs)
{
    forAll(ibs, ibi)
    {
        ibs[ibi].getGeomModel().getBBoxes();
    }

    // use insertion sort to sort the neighbour list
    // besides the first sorting, this is O(N) efficient
    for (label coord = 0; coord < 3; ++coord)
    {
        if (!verletLists_[coord].empty())
        {
            auto it2 = verletLists_[coord].begin();
            auto it1 = it2++;

            while(true)
            {
                if((*it1)->getPoint()[coord]
                    > (*it2)->getPoint()[coord])
                {
                    swapVerletPoints(*it1, *it2, coord);
                    std::swap(*it1, *it2);

                    if (it1 != verletLists_[coord].begin())
                    {
                        it2 = it1--;
                        continue;
                    }
                }
                it1 = it2++;
                if (it2 == verletLists_[coord].end())
                {
                    break;
                }
            }
        }
    }
}
//---------------------------------------------------------------------------//
