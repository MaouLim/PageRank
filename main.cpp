#include <iostream>
#include <fstream>

#include <spare_matrix_f.h>
#include <solver.h>

static const std::string data_file = "graph.data";

void handle_err(const std::string& msg) {
	std::cout << "FATAL: " << msg << std::endl;
	std::cout << "Enter a key to quit..." << std::endl;
	std::cin.get();
	exit(-1);
}

std::string get_url(const std::string& url_file, size_t id) {
	std::ifstream stream(url_file, std::ios::in);
	assert(stream.is_open());

	std::string line;

	size_t i = 0;
	while (i <= id) { std::getline(stream, line); ++i; }

	std::stringstream sstream(line);

	size_t seq = 0;
	sstream >> seq;
	assert(seq == id);

	std::string url;
	sstream >> url;

	return url;
}

int main(int argc, char** argv) {
	if (3 != argc) { handle_err("Invalid console parameters"); }

	const std::string url_file = argv[1];
	const std::string out_file = argv[2];

	std::cout << "Start loading matrix..." << std::endl;

	if (!tools::convert<size_t, double>(url_file, data_file)) {
		handle_err("Failed to load matrix");
	}

	tools::spare_matrix_f<size_t, double> mat(data_file);

	std::cout << "Finish loading matrix." << std::endl;

	std::cout << "Start iterating..." << std::endl;

	page_rank::page_rank_solver_f solver(0.2, 200, 1e-8);
	auto result = solver.solve(mat);

	std::cout << "Finish iterating." << std::endl;

	std::cout << "Iteration: " << result.first << std::endl;

	/* convert to probability distribution */
	const auto sum = tools::sum(result.second);
	for (auto& each : result.second) { each /= sum; }

	std::cout << "Start sorting..." << std::endl;

	std::vector<size_t> indices(result.second.size(), 0);
	for (size_t i = 0; i < indices.size(); ++i) { indices[i] = i; }

	std::sort(
		indices.begin(), indices.end(),
		[&result](size_t a, size_t b) { return result.second[a] > result.second[b]; }
	);

	std::cout << "Finish sorting." << std::endl;

	std::ofstream stream(out_file, std::ios::out);
	if (!stream.is_open()) { handle_err("Failed to open output file"); }

	for (size_t i = 0; i < 10; ++i) {
		stream << get_url(url_file, indices[i]) << ' ' 
		       << result.second[indices[i]]     << std::endl;
	}

	std::cout << "Task finished, enter a key to quit..." << std::endl;
	std::cin.get();
	return 0;
}