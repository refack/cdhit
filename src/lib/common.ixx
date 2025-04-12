export module common;

import std;

// inline some macros
#ifdef _CDHIT_VERSION_
export constexpr std::string_view CDHIT_VERSION{_CDHIT_VERSION_};
#else
export constexpr std::string_view CDHIT_VERSION{"0.0.0"};
#endif
#undef _CDHIT_VERSION_

#ifdef OPENMP
export constexpr std::string_view WITH_OPENMP = "(+OpenMP)";
#else
export constexpr std::string_view WITH_OPENMP{};
#endif
// end inlining

export constexpr std::string_view BUILD_DATE{__DATE__};

export constexpr size_t MAX_UAA = 4;
export constexpr size_t MAX_FILE_NAME = 1280;
export constexpr size_t MAX_SEG = 50;
export constexpr size_t MAX_BIN_SWAP = 2 << 30;
export constexpr size_t MAX_TABLE_SIZE = 5'0000'000;
export constexpr size_t CLOCK_TICKS = 100;
export constexpr size_t FAILED_FUNC = 1;
export constexpr size_t OK_FUNC = 0;

export consteval void SetGlobals(const size_t uaa) {
    // MAX_UAA = uaa;
}




#if defined(LONG_SEQ) && defined(SHORT_SEQ)
#error "LONG_SEQ and SHORT_SEQ cannot be defined at the same time"
#endif

// class function definition
// constexpr char aa[] = { "ARNDCQEGHILKMFPSTWYVBZX" };
// {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,2,6,20};

export constexpr size_t MEGA_MiBi = 1'000'000;

#if defined(LONG_SEQ)
// if the longset sequence is longer than 65535, use uint32_t
export typedef std::uint32_t INTs;
#elif defined(SHORT_SEQ)
export typedef std::uint8_t INTs;
#else
export typedef std::uint16_t INTs;
#endif
export typedef std::vector<INTs> VectorIntX;
export typedef std::vector<VectorIntX> MatrixIntX;

export typedef std::vector<size_t> VectorInt64;
export typedef VectorInt64 VectorInt;
export typedef VectorInt64 VectorInt32;

export typedef std::vector<VectorInt64> MatrixInt64;
export typedef MatrixInt64 MatrixInt;

export constexpr size_t NAA0 = 1;
export constexpr size_t NAA1 = MAX_UAA;
export constexpr size_t NAA2 = NAA1 * NAA1;
export constexpr size_t NAA3 = NAA1 * NAA2;
export constexpr size_t NAA4 = NAA2 * NAA2;
export constexpr size_t NAA5 = NAA2 * NAA3;
export constexpr size_t NAA6 = NAA3 * NAA3;
export constexpr size_t NAA7 = NAA3 * NAA4;
export constexpr size_t NAA8 = NAA4 * NAA4;
export constexpr size_t NAA9 = NAA4 * NAA5;
export constexpr size_t NAA10 = NAA5 * NAA5;
export constexpr size_t NAA11 = NAA5 * NAA6;
export constexpr size_t NAA12 = NAA6 * NAA6;
export constexpr std::array NAAN_array{
	NAA0,
	NAA1,
	NAA2,
	NAA3,
	NAA4,
	NAA5,
	NAA6,
	NAA7,
	NAA8,
	NAA9,
	NAA10,
	NAA11,
	NAA12,
};

// either 4'000'000 or (1<<22)
export constexpr auto MAX_TABLE_SEQ = 4'000'000;

export enum
{
	DP_BACK_NONE = 0,
	DP_BACK_LEFT_TOP = 1,
	DP_BACK_LEFT = 2,
	DP_BACK_TOP = 3
};

export struct IndexCount
{
    size_t index = 0;
    size_t count = 0;
};

export typedef std::vector<IndexCount> NVectorIndexCount;

export template<class T>
class VectorIterFactory {
    std::vector<T>& vec;
    const size_t frag_size;

public:
    explicit VectorIterFactory(std::vector<T>& v, const size_t _frag_size) : vec(v), frag_size(_frag_size) {}

    constexpr typename std::vector<T>::iterator operator()(const std::ptrdiff_t fi) const {
        auto i = fi * frag_size;
        if (i >= vec.size()) return vec.end();
        return vec.begin() + i;
    }
};


