/*!
 * @file overlap.hpp
 *
 * @brief Overlap class header file
 */

#include "pile.hpp"
#include "overlap.hpp"

namespace rala {

Overlap::Overlap(uint64_t a_id, uint64_t b_id, double, uint32_t,
    uint32_t a_rc, uint32_t a_begin, uint32_t a_end, uint32_t a_length,
    uint32_t b_rc, uint32_t b_begin, uint32_t b_end, uint32_t b_length)
        : a_name_(), a_id_(a_id - 1), a_begin_(a_begin), a_end_(a_end),
        a_length_(a_length), b_name_(), b_id_(b_id - 1), b_begin_(b_begin),
        b_end_(b_end), b_length_(b_length), length_(std::max(a_end - a_begin,
        b_end - b_begin)), orientation_(a_rc == b_rc ? 0 : 1),
        is_transmuted_(false) {
}

Overlap::Overlap(const char* a_name, uint32_t a_name_length,
    uint32_t a_length, uint32_t a_begin, uint32_t a_end, char orientation,
    const char* b_name, uint32_t b_name_length, uint32_t b_length,
    uint32_t b_begin, uint32_t b_end, uint32_t, uint32_t overlap_length, uint32_t)
        : a_name_(a_name, a_name_length), a_id_(), a_begin_(a_begin),
        a_end_(a_end), a_length_(a_length), b_name_(b_name, b_name_length),
        b_id_(), b_begin_(b_begin), b_end_(b_end), b_length_(b_length),
        length_(overlap_length), orientation_(orientation == '+' ? 0 : 1),
        is_transmuted_(false) {
}

Overlap::~Overlap() {
}

bool Overlap::transmute(const std::vector<std::unique_ptr<Pile>>& piles,
    const std::unordered_map<std::string, uint64_t>& name_to_id) {

    if (is_transmuted_) {
        return true;
    }

    if (!a_name_.empty()) {
        auto it = name_to_id.find(a_name_);
        if (it == name_to_id.end()) {
            return false;
        }
        a_id_ = it->second;
        std::string().swap(a_name_);
    }
    if (a_id_ >= piles.size() || piles[a_id_] == nullptr) {
        return false;
    }
    if (a_length_ != piles[a_id_]->data().size()) {
        fprintf(stderr, "[rala::Overlap::transmute] error: "
            "unequal lengths in sequence and overlap file for sequence with id %lu!\n",
            a_id_);
        exit(1);
    }

    if (!b_name_.empty()) {
        auto it = name_to_id.find(b_name_);
        if (it == name_to_id.end()) {
            return false;
        }
        b_id_ = it->second;
        std::string().swap(b_name_);
    }

    if (b_id_ >= piles.size() || piles[b_id_] == nullptr) {
        return false;
    }
    if (b_length_ != piles[b_id_]->data().size()) {
        fprintf(stderr, "[rala::Overlap::transmute] error: "
            "unequal lengths in sequence and overlap file for sequence with id %lu!\n",
            b_id_);
        exit(1);
    }

    is_transmuted_ = true;
    return true;
}

bool Overlap::trim(const std::vector<std::unique_ptr<Pile>>& piles) {

    if (!is_transmuted_) {
        fprintf(stderr, "[racon::Overlap::trim] error: overlap is not transmuted!\n");
        exit(1);
    }
    if (a_id_ >= piles.size() || piles[a_id_] == nullptr ||
        b_id_ >= piles.size() || piles[b_id_] == nullptr) {
        return false;
    }

    const auto& pile_a = piles[a_id_];
    const auto& pile_b = piles[b_id_];

    if (pile_a->begin() > a_length_ || pile_a->end() > a_length_ ||
        pile_b->begin() > b_length_ || pile_b->end() > b_length_) {

        fprintf(stderr, "[rala::Overlap::trim] error: "
            "invalid trimmed begin, end coordinates!\n");
        exit(1);
    }

    if (a_begin_ >= pile_a->end() || a_end_ <= pile_a->begin()) {
        return false;
    }
    if (b_begin_ >= pile_b->end() || b_end_ <= pile_b->begin()) {
        return false;
    }

    uint32_t a_new_begin = a_begin_ + (b_begin_ < pile_b->begin() ?
        pile_b->begin() - b_begin_ : 0);
    uint32_t a_new_end = a_end_ - (b_end_ > pile_b->end() ?
        b_end_ - pile_b->end() : 0);
    if (a_new_begin >= pile_a->end() || a_new_end <= pile_a->begin()) {
        return false;
    }

    uint32_t b_new_begin = b_begin_ + (a_begin_ < pile_a->begin() ?
        pile_a->begin() - a_begin_ : 0);
    uint32_t b_new_end = b_end_ - (a_end_ > pile_a->end() ?
        a_end_ - pile_a->end() : 0);
    if (b_new_begin >= pile_b->end() || b_new_end <= pile_b->begin()) {
        return false;
    }

    a_new_begin = std::max(a_new_begin, pile_a->begin()) - pile_a->begin();
    a_new_end = std::min(a_new_end, pile_a->end()) - pile_a->begin();
    b_new_begin = std::max(b_new_begin, pile_b->begin()) - pile_b->begin();
    b_new_end = std::min(b_new_end, pile_b->end()) - pile_b->begin();

    if (a_new_begin >= a_new_end || b_new_begin >= b_new_end) {
        return false;
    }

    a_begin_ = a_new_begin;
    a_end_ = a_new_end;
    a_length_ = pile_a->end() - pile_a->begin();

    b_begin_ = b_new_begin;
    b_end_ = b_new_end;
    b_length_ = pile_b->end() - pile_b->begin();

    length_ = std::max(a_end_ - a_begin_, b_end_ - b_begin_);

    return true;
}

OverlapType Overlap::type() const {

    uint32_t a_begin = a_begin_;
    uint32_t a_end = a_end_;
    uint32_t b_begin = orientation_ == 0 ? b_begin_ : b_length_ - b_end_;
    uint32_t b_end = orientation_ == 0 ? b_end_ : b_length_ - b_begin_;

    uint32_t overhang = std::min(a_begin, b_begin) + std::min(a_length_ -
        a_end, b_length_ - b_end);

    if (a_end - a_begin < (a_end - a_begin + overhang) * 0.875 ||
        b_end - b_begin < (b_end - b_begin + overhang) * 0.875) {
        return OverlapType::kX;
    }
    if (a_begin <= b_begin && (a_length_ - a_end) <= (b_length_ - b_end)) {
        return OverlapType::kB;
    }
    if (a_begin >= b_begin && (a_length_ - a_end) >= (b_length_ - b_end)) {
        return OverlapType::kA;
    }

    auto absolute_difference = [](uint32_t a, uint32_t b) -> uint32_t {
        return a > b ? (a - b) : (b - a);
    };

    if (absolute_difference(a_end_ - a_begin_, b_end_ - b_begin_) < length_ * 0.01) {
        uint32_t min_extension = 0.05 * std::max(a_length_, b_length_);

        if (absolute_difference(a_begin, b_begin) < min_extension) {
            if ((a_length_ - a_end) >= (b_length_ - b_end)) {
                return OverlapType::kA;
            } else {
                return OverlapType::kB;
            }
        }
        if (absolute_difference((a_length_ - a_end), (b_length_ - b_end)) < min_extension) {
            if (a_begin >= b_begin) {
                return OverlapType::kA;
            } else {
                return OverlapType::kB;
            }
        }
    }

    if (a_begin > b_begin) {
        return OverlapType::kAB;
    }
    return OverlapType::kBA;
}

}
