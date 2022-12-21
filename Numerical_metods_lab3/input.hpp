#pragma once
#include <cinttypes>
#include <json.hpp>


struct input_t
{
	double eps;
	uint32_t n, m;
	uint64_t max_step;
};

void get_input(const nlohmann::json& from, input_t& input);