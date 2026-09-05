#include <sys/ioctl.h>
#include <unistd.h>
#include <chrono>
#include <thread>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>



class progressBar {
public:
    progressBar(volatile double* progress, const int nProgress, const int width, const std::string* title) {
        constexpr int offset = 19;

        _progress = progress;
        _nProgress = nProgress;

        if (title) {
            _title = *title + ": ";
        }

        int _width;
        if (width) {
            _width = width;
        }
        else {
            winsize win{};
            ioctl(STDOUT_FILENO, TIOCGWINSZ, &win);
            _width = win.ws_col;
        }
        _width = std::max(_width, 80);

        _barLength = _width - static_cast<int>(_title.length()) - offset;

        _start = std::chrono::steady_clock::now();
    }

    ~progressBar() {
        std::cout << std::endl;
    }

    void show() const {
        double minProgress = *std::ranges::min_element(_progress, _progress + _nProgress);
        minProgress = std::min(1.0, std::max(static_cast<double>(minProgress), 0.0));
        int minBarProgress = static_cast<int>(_barLength * minProgress);

        double avgProgress = 0;
        for (int i = 0; i < _nProgress; i++) {
            avgProgress += _progress[i] / _nProgress;
        }
        avgProgress = std::min(1.0, std::max(static_cast<double>(avgProgress), 0.0));
        int avgBarProgress = static_cast<int>(_barLength * avgProgress);
        double avgBarFraction = _barLength * avgProgress - avgBarProgress;

        std::string bar;
        for (int i = 0; i < minBarProgress; i++) {
            bar += "\u2591";
        }
        for (int i = 0; i < avgBarProgress - minBarProgress; i++) {
            bar += "\u2588";
        }
        if (avgBarFraction > 0.0) {
            bar += _barChars[static_cast<size_t>(avgBarFraction * 8)];
            avgBarProgress++;
        }
        for (int i = 0; i < _barLength - avgBarProgress; i++) {
            bar += ' ';
        }

        const auto elapsed = duration_cast<std::chrono::seconds>((std::chrono::steady_clock::now() - _start));
        const std::chrono::hh_mm_ss time{elapsed};

        std::cout << "\r" << _title
                  << std::setw(5) << std::fixed << std::setprecision(1) << (avgProgress * 100.0)
                  << "%|" << bar << "| [" << time << "]" << std::flush;
    }

private:
    static constexpr std::string _barChars[8] {
        " ",      "\u258F",
        "\u258E", "\u258D",
        "\u258C", "\u258B",
        "\u258A", "\u2589"
    };

    volatile double* _progress;
    int _nProgress;

    std::string _title{};

    int _barLength;

    std::chrono::time_point<std::chrono::steady_clock> _start;
};



int main() {
    volatile double* progress = new double[16]();

    constexpr std::string title = "Progress";
    const auto pb = progressBar(progress, 16, 0, &title);

    std::vector<std::jthread> threads;
    threads.reserve(16);

    for (int i = 0; i < 16; i++) {
        threads.emplace_back([progress, i]() {
            while (progress[i] < 1.0) {
                progress[i] += (std::rand() % 100 - 40) / 10000.0;
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
            }
        });
    }

    while (*std::ranges::min_element(progress, progress + 16) < 1.0) {
        pb.show();

        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
    pb.show();

    delete[] progress;
    return 0;
}
