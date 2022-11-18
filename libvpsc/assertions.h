#pragma once
#define COLA_ASSERT(expr)  static_cast<void>(0)
#include <cstring>
#include <sstream>
#include <cassert>
namespace vpsc 
{
	class CriticalFailure
	{
	public:
		CriticalFailure(const char* expr, const char* file, int line,
			const char* function = nullptr)
			: expr(expr),
			file(file),
			line(line),
			function(function)
		{
		}
		std::string what() const
		{
			std::stringstream s;
			s << "ERROR: Critical assertion failed.\n";
			s << "  expression: " << expr << "\n";
			s << "  at line " << line << " of " << file << "\n";
			if (function)
			{
				s << "  in: " << function << "\n";
			}

			return s.str();
		}
	private:
		const char* expr;
		const char* file;
		int line;
		const char* function;
	};

}


