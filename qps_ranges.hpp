/*
 *       /\
 *      /__\       Matteo Cicuttin (C) 2016 - matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code for scientific publications, you are required to
 * cite it.
 */

#pragma once

namespace disk {

class qps_range
{
    size_t m_min, m_max;

public:
    qps_range()
        : m_min(0), m_max(0)
    {}

    qps_range(size_t min, size_t max)
        : m_min(min), m_max(max)
    {
        assert(m_min <= m_max);
    }

    size_t from() const {
        return m_min;
    }

    size_t min() const {
        return m_min;
    }

    size_t to() const {
        return m_max;
    }

    size_t max() const {
        return m_max;
    }

    size_t size() const {
        return m_max - m_min;
    }

    qps_range remove_offset() const {
        return qps_range(0, m_max - m_min);
    }

    bool operator<(const qps_range& other) // tell if this range is contained in other
    {
        return (m_min >= other.m_min) && (m_max <= other.m_max);
    }
};

std::ostream&
operator<<(std::ostream& os, const qps_range& range)
{
    os << "[" << range.min() << ", " << range.max() << ") " << range.size();
    return os;
}

class qpspace_ranges
{
    qps_range                   m_cell_range;
    qps_range                   m_face_range;
    size_t                      m_num_faces;

public:
    qpspace_ranges()
        : m_num_faces(0)
    {}

    qpspace_ranges(size_t num_cell_qps, size_t num_face_qps, size_t num_faces)
    {
        m_cell_range = qps_range(0, num_cell_qps);
        m_face_range = qps_range(0, num_face_qps);
        m_num_faces = num_faces;
    }

    qps_range cell_range(void) const {
        return m_cell_range;
    }

    qps_range face_range(size_t face) const {
        assert(face < m_num_faces);
        size_t face_range_size  = m_face_range.size();
        return qps_range( face*face_range_size, (face+1)*face_range_size);
    }

    qps_range all_faces_range(void) const {
        size_t face_start = 0;
        size_t face_end   = m_num_faces * m_face_range.size();
        return qps_range(face_start, face_end);
    }

    size_t total_size(void) const {
        return m_cell_range.size() + m_num_faces * m_face_range.size();
    }

    size_t num_faces(void) const {
        return m_num_faces;
    }
};

std::ostream&
operator<<(std::ostream& os, const qpspace_ranges& dsr)
{
    os << "Cell range : " << dsr.cell_range() << std::endl;
    os << "Face ranges: " << std::endl;
    for (size_t i = 0; i < dsr.num_faces(); i++)
        os << " - " << dsr.face_range(i) << std::endl;

    return os;
}


template<typename T>
std::vector<T>
take(const std::vector<T>& from, const qps_range& range)
{
    auto begin  = std::advance(from.begin(), range.from());
    auto end    = std::advance(from.end(), range.to());
    return std::vector<T>(begin, end);
}

template<typename T>
dynamic_vector<T>
take(const dynamic_vector<T>& from, const qps_range& range)
{
    return from.head(range.to()).tail(range.size());
}

template<typename T>
/*dynamic_matrix<T>*/auto
take(const dynamic_matrix<T>& from, const qps_range& row_range,
     const qps_range& col_range)
{
    return from.block(row_range.from(), col_range.from(),
                      row_range.size(), col_range.size());
}

} //namespace disk
