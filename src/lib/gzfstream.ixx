module;

#include <zlib.h>

export module gzfstream;

import std;

namespace HitFileStream {
using std::make_unique;

export constexpr size_t MAX_LINE_SIZE = 300'000;
export constexpr size_t LINE_BUF_SIZE = MAX_LINE_SIZE + 1;


class Wrapper {
public:
    virtual ~Wrapper() = default;
    virtual void open(const char* filename) = 0;
    virtual void close() = 0;
    virtual size_t read(char* buffer, std::streamsize size) = 0;
    virtual void seek(z_off_t pos) = 0;
};


class CFile_Wrapper : public Wrapper {
private:
    std::ifstream file_;

public:
    void open(const char* filename) override {
        file_.open(filename, std::ios::binary);
        if (!file_)
            throw std::invalid_argument("Failed to open input file");
    }

    void close() override { file_.close(); }

    size_t read(char* buffer, std::streamsize size) override {
        file_.read(buffer, size);
        return file_.gcount();
    }

    void seek(z_off_t pos) override {
        file_.seekg(pos);
        if (!file_)
            throw std::runtime_error("seekg failed");
    }
};


class GZFile_Wrapper : public Wrapper {
public:
    void open(const char* filename) override {
        file_.reset(gzopen(filename, "rb"));
        if (!file_)
            throw std::invalid_argument("Failed to open gzipped input file");
    }

    void close() override {
        if (file_)
            gzclose(file_.release());
    }

    size_t read(char* read_buffer, std::streamsize read_size) override { return gzread(file_.get(), read_buffer, read_size); }

    void seek(z_off_t pos) override {
        if (gzseek(file_.get(), pos, SEEK_SET) < 0)
            throw std::runtime_error("gzseek failed");
    }

private:
    struct gzCloser {
        void operator()(gzFile f) const {
            if (f)
                gzclose(f);
        }
    };

    std::unique_ptr<gzFile_s, gzCloser> file_;
};


class FileWrapper_StreamBuffer : public std::streambuf {
private:
    Wrapper* file_;
    std::array<char, MAX_LINE_SIZE> buffer_{};

public:
    explicit FileWrapper_StreamBuffer(Wrapper* file) : file_(file) { setg(buffer_.data(), buffer_.data(), buffer_.data()); }

protected:
    int_type underflow() override {
        if (!file_)
            return traits_type::eof();
        size_t bytes_read = file_->read(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
        if (bytes_read <= 0)
            return traits_type::eof();
        setg(buffer_.data(), buffer_.data(), buffer_.data() + bytes_read);
        return traits_type::to_int_type(*gptr());
    }
};


export class HitFStream : public std::istream {
private:
    std::unique_ptr<Wrapper> file_;
    std::unique_ptr<FileWrapper_StreamBuffer> buf_;

public:
    explicit HitFStream(const std::string_view filename) : std::istream(nullptr), file_(nullptr), buf_(nullptr) {
        if (filename.ends_with(".gz"))
            file_ = make_unique<GZFile_Wrapper>();
        else
            file_ = make_unique<CFile_Wrapper>();
        file_->open(filename.data());
        buf_ = make_unique<FileWrapper_StreamBuffer>(file_.get());
        rdbuf(buf_.get());
    }

    void seek(z_off_t pos) {
        file_->seek(pos);
        buf_ = std::make_unique<FileWrapper_StreamBuffer>(file_.get());
        rdbuf(buf_.get());
    }
};
} // namespace HitFileStream
