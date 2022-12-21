#include "input.hpp"


using nlohmann::json;

void from_json(const json& j, input_t& input)
{
	input.n = j.at("n");
	input.m = j.at("m");

	input.max_step = j.at("max_step");
	input.eps = j.at("eps");
}



void get_input(const json& from, input_t& input)
{
	input = from.at("input");
}