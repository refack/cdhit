export module Sequence;

import std;

import common;
import Options;


export enum class seq_state {
	NONE = 0,
	IS_REP = 1,
	IS_REDUNDANT = 2,
	IS_PROCESSED = 16,
	IS_MINUS_STRAND = 32
};

export constexpr bool operator&(const seq_state& lhs, const seq_state& rhs) {
	return 
		static_cast<std::underlying_type_t<seq_state>>(lhs) &
		static_cast<std::underlying_type_t<seq_state>>(rhs)
		;
}

export constexpr seq_state& operator|=(seq_state& lhs, seq_state rhs) {
	lhs = static_cast<seq_state>(
		static_cast<std::underlying_type_t<seq_state>>(lhs) |
		static_cast<std::underlying_type_t<seq_state>>(rhs)
	);
	return lhs;
}

export constexpr seq_state& operator-=(seq_state& lhs, seq_state rhs) {
	lhs = static_cast<seq_state>(
		static_cast<std::underlying_type_t<seq_state>>(lhs) &
		~static_cast<std::underlying_type_t<seq_state>>(rhs)
	);
	return lhs;
}


// static_assert(std::is_move_constructible_v<Sequence>);
// static_assert(std::is_move_assignable_v<Sequence>);
export struct Sequence {
private:
    // real sequence, if it is not stored swap file:
    std::vector<char> _data{};
    std::string _identifier{};


    // uint32_t stats;

    // // if swap != nullptr, the sequence is stored in file.
    // // swap is opened as temporary file, which will be deleted automatically
    // // after the program is finished:
    // std::ifstream swap{};
    // // stream offset of the sequence:
    // size_t offset = 0;
public:
    // length of the sequence:
    size_t size = 0;
    // allocated size of the data buffer:
    size_t bufsize = 0;

    // stream offset of the description string in the database (for swap file):
    size_t des_begin = 0;
    // total record length (for swap file):
    size_t tot_length = 0;

    // when we store two back-to-back merged sequences:
    size_t size_R2 = 0; // for back-to-back merged seq size.R1 = size - size.R2
    size_t des_begin2 = 0;
    size_t tot_length2 = 0;

    // index of the sequence in the original database:
    size_t index = 0;
    size_t cluster_id = 0;
    float identity = 0.0f;
    float distance = 2.0f;
    size_t coverage[4] = {0, 0, 0, 0};
    seq_state state{};

    [[nodiscard]] size_t get_size_R1() const { return size - size_R2; }

    // #####################################
    // ########### temporary API ###########
    // #####################################

    // Getter
    char* _get_data() { return _data.data(); }
    // Property declaration (MSVC extension)
    __declspec(property(get = _get_data)) char* data;

    // Getter
    char* _get_identifier() { return _identifier.data(); }
    // Property declaration (MSVC extension)
    __declspec(property(get = _get_identifier)) char* identifier;

    // #####################################
    // ######### end temporary API #########
    // #####################################

    // Non-copyable since no use case for that
    Sequence(Sequence& other) = delete;
    Sequence& operator=(const Sequence& other) = delete;

    Sequence() = default;

    Sequence(const Sequence& other) = default;

    Sequence(Sequence&& other) = default;

    Sequence& operator=(Sequence&& other) = default;

    // Should only be used when loading PE from database
    Sequence(Sequence&& R1, Sequence&& R2): Sequence(std::move(R1)) {
        // back to back merge for pair-ends
        //
        // PE sequenced fragments
        // R1 -> XXXXXXABC ------------------- NMLYYYYYY <--R2
        //
        // from the database
        // >R1           >R2
        // XXXXXXABC     YYYYYYLMN
        //
        // what we want
        // >R12
        // NMLYYYYYYXXXXXXABC

        // update info for R2
        size += R2.size;
        bufsize += R2.size;
        size_R2 = R2.size;

        des_begin2 = R2.des_begin;
        tot_length2 = R2.tot_length;

        // The rest of the params should still be 0

        // hold the R1 seq
        auto tmp = _data;
        // move the revered R2 first
        _data = std::move(R2._data);
        // inplace reverse R2 in _data
        std::ranges::reverse(_data);
        // concat forward R1
        _data.insert_range(_data.end(), tmp);
    }

    ~Sequence() = default;

    [[nodiscard]] auto get_data() { return std::string_view{_data}; }
    [[nodiscard]] const auto get_data() const { return std::string_view{_data}; }
    [[nodiscard]] auto get_identifier() const { return std::string_view{_identifier}; }
    [[nodiscard]] size_t get_data_size() const { return _data.size(); }
    [[nodiscard]] bool id_empty() const { return _identifier.empty(); }

    void Clear() {
        _data.clear();
        _identifier.clear();
        bufsize = 0;
        size = 0;
    }

    void Resize(size_t new_capacity) {
        if (new_capacity == bufsize)
            return;
        _data.resize(new_capacity + 1, '\0');
        bufsize = _data.size();
        size = new_capacity;
    }

    void Reserve(size_t new_size) {
        size = new_size;
        if (size > bufsize) {
            _data.reserve(size + 1);
            _data.back() = '\0';
            bufsize = _data.size();
        }
    }

    Sequence& operator=(const char* s) {
        size = 0; // avoid copying;
        size_t other_len = std::strlen(s);
        Resize(other_len);
        std::memcpy(data, s, other_len + 1);
        return *this;
    }

    bool operator<(const Sequence& other) const { return this->_data.size() < other._data.size(); }

    void operator+=(const char* s) {
        size_t m = size;
        size_t n = std::strlen(s);
        Reserve(m + n);
        std::memcpy(data + m, s, n);
    }

    void operator+=(const std::string& s) {
        size_t m = size;
        size_t n = s.size();
        Reserve(m + n);
        std::memcpy(data + m, s.c_str(), n);
    }

    [[nodiscard]] size_t cutoff2_len() const { return std::llround(static_cast<double>(size) * options.diff_cutoff2); }

    void set_id(Sequence& other) {
        this->_identifier = std::move(other._identifier);
        other._identifier.clear();
    }
    void set_id(std::string&& other) { _identifier = std::move(other); }

    bool FormatAndValidate() // FormatAndValidate
    {
        for (char& ch : std::views::reverse(_data)) {
            if (std::isspace(ch) || ch == '*') {
                size--;
                ch = '\0';
            }
        }
        _data.resize(size);
        for (char& ch : _data) {
            bool is_valid = std::isalpha(ch) || std::isspace(ch);
            if (!is_valid)
                return true;
            ch = static_cast<char>(std::toupper(ch));
        }
        return id_empty();
    }

    void ConvertBases() {
        for (size_t i = 0; i < size; i++)
            data[i] = aa2idx[data[i] - 'A'];
    }


    void trim(size_t trim_len) {
        if (trim_len >= size)
            return;
        size = trim_len;
        if (size)
            data[size] = '\0';
    }


    friend std::ostream& operator<<(std::ostream& fout, const Sequence& seq) {
        auto bare_id = std::string_view{seq._identifier.data() + 1, seq._identifier.size() - 1};
        std::print(fout, "{}\t{}{}, >{}...", seq.cluster_id, seq.size, (options.isEST ? "nt" : "aa"), bare_id);

        if (seq.id_empty())
            return fout << " *\n";

        fout << " at ";
        if (options.print != 0) {
            std::print(fout, "{}:{}:{}:{}/", seq.coverage[0], seq.coverage[1], seq.coverage[2], seq.coverage[3]);
        }
        if (options.isEST) {
            std::print(fout, "{}/", (seq.state & seq_state::IS_MINUS_STRAND) ? '-' : '+');
        }
        std::print(fout, "{:.2f}%", seq.identity * 100);
        if (options.useDistance()) {
            std::print(fout, "/{:.2f}%", seq.distance * 100);
        }

        return fout << std::endl;
    }
};
static_assert(std::is_swappable_v<Sequence>);
static_assert(std::indirectly_writable<std::vector<Sequence>::iterator, Sequence>);


// void Sequence::SwapIn()
// {
// 	std::ifstream ifs;
// 	auto filebuf = ifs.rdbuf();
// 	std::cref(*filebuf);
// 	ifs.native_handle();
// 	if (data) return;
// 	if (swap == nullptr) throw std::runtime_error("Can not swap in sequence");
// 	Resize(size);
// 	std::fseek(swap, offset, std::SEEK_SET);
// 	if (std::fread(data, 1, size, swap) == 0) throw std::runtime_error("Can not swap in sequence");
// 	data[size] = '\0';
// }
//
// void Sequence::SwapOut()
// {
// 	if (swap && _data) {
// 		_data.reset();
// 		bufsize = 0;
// 	}
// }
