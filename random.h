#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include <chrono>

namespace Random {
	// Returns a seeded Mersenne Twister
	// Note: we'd prefer to return a std::seed_seq (to initialize a std::mt19937), but std::seed can't be copied, so it can't be returned by value.
	// Instead, we'll create a std::mt19937, seed it, and then return the std::mt19937 (which can be copied).
	inline std::mt19937 generate() {
		std::random_device rd{};

		// Create seed_seq with clock and 7 random numbers from std::random_device
		std::seed_seq ss{
			static_cast<std::seed_seq::result_type>(std::chrono::steady_clock::now().time_since_epoch().count()),
				rd(), rd(), rd(), rd(), rd(), rd(), rd() };

		return std::mt19937{ ss };
	}

	// Here's our global std::mt19937 object.
	// The inline keyword means we only have one global instance for our whole program.
	thread_local inline std::mt19937 mt{ generate() }; // generates a seeded std::mt19937 and copies it into our global object

	// Generate a random double between [min, max] (inclusive)
	inline double uniform_double(double min, double max) {
		return std::uniform_real_distribution<double>{min, max}(mt);
	}

    // Generate a random double X ~ N(mean, stddev)
	inline double gaussian_double(double mean, double stddev) {
	    return std::normal_distribution<double>{mean, stddev}(mt);
	}

    // Generate a random complex X + iY with X,Y ~ N(mean, stddev)
	inline std::complex<double> gaussian_complex(double mean, double stddev) {
		thread_local std::normal_distribution<double> dist{mean, stddev};
	    double x = dist(mt);
	    double y = dist(mt);

	    return {x,y};
	}
}

#endif //RANDOM_H
