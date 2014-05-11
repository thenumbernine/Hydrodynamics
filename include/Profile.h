#pragma once

#include "Common/Exception.h"

#include <math.h>
#include <sys/time.h>

#include <sstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <functional>
#include <algorithm>

struct Profiler {
	static double getTime() {
		timeval t;
		gettimeofday(&t, NULL);
		return (double)t.tv_sec + (double)t.tv_usec * 1e-6;
	}

	//taken from universe/offline/Stat
	struct Stat {
		double v, n;
		double avg, sqavg;
		double min, max;
		Stat() : v(0), n(0), avg(0), sqavg(0), min(HUGE_VAL), max(-HUGE_VAL) {}
		void accum(double v) {
			++n;
			min = std::min(min, v);
			max = std::max(max, v);
			avg += (v - avg) / n;
			sqavg += (v * v - sqavg) / n;
		}
		double stddev() { 
			return sqrt(sqavg - avg * avg);
		}
		std::string str() {
			std::ostringstream o;
			o << "[" << min << " ... " << avg << " +- " << stddev() << " ... " << max << "]";
			return o.str();
		}
	};

	struct ProfileEntry {
		Stat duration;
		Stat durationExcludingChildren;
		const char *name;
		const char *parentName;
		int order;

		ProfileEntry()
		: name(NULL)
		, parentName(NULL)
		, order(0)
		{}

		std::map<const char *, bool> childNames;

		void print(int spaces = 0) {
			std::cout << "#" << order << "\t";
			for (int i = 0; i < spaces; ++i) std::cout << " ";
			std::cout 
				<< "duration " << durationExcludingChildren.str()
				<< " inclusive " << duration.str()
				<< " \t"
				<< name << std::endl;
			++spaces;	
			std::for_each(childNames.begin(), childNames.end(), [&](std::map<const char *, bool>::value_type v) {
				profileMap[v.first].print(spaces);
			});
		}
	};

	static std::map<const char *, ProfileEntry> profileMap;
	
	const char *name;
	double startTime, subDurations;

	Profiler *lastProfiler;
	static Profiler *currentProfiler;
	static const char *rootName;	//cache for done()
	
	Profiler(const char *name_) : name(name_) {
		lastProfiler = currentProfiler;
		currentProfiler = this;
		startTime = getTime();
		subDurations = 0;
	}

	~Profiler() {
		double endTime = getTime();
		double duration = endTime - startTime;
		if (lastProfiler) {
			lastProfiler->subDurations += duration;
		}

		rootName = name;
		currentProfiler = lastProfiler;

		ProfileEntry &entry = profileMap[name];
		entry.name = name;
		entry.duration.accum(duration);
		entry.durationExcludingChildren.accum(duration - subDurations);
		if (lastProfiler) {
			entry.parentName = lastProfiler->name;
			profileMap[lastProfiler->name].childNames[entry.name] = true;
		}
	}

	static void done() {
		if (!rootName) throw Exception() << "need to run Profiler at least once";
		
		//sort profiles by minimum duration excluding children
		std::vector<std::pair<const char *, double> > profileOrder;
		std::for_each(profileMap.begin(), profileMap.end(), [&](const std::pair<const char *, ProfileEntry> &pair) {
			profileOrder.push_back(std::pair<const char *, double>(pair.first, pair.second.durationExcludingChildren.min));
		});
		std::sort(profileOrder.begin(), profileOrder.end(), [&](std::pair<const char *, double> a, std::pair<const char *, double> b) {
			return a.second > b.second;
		});
		int order = 0;
		std::for_each(profileOrder.begin(), profileOrder.end(), [&](std::pair<const char *, double> pair){
			profileMap[pair.first].order = ++order;
		});
		
		profileMap[rootName].print();
	}
};


//http://stackoverflow.com/questions/3859944/combining-string-literals-and-integer-constants
#define STRINGIZE_DETAIL_(v)	#v
#define STRINGIZE(v) STRINGIZE_DETAIL_(v)
#define FUNCTION_STRING	STRINGIZE(__func__)	//not working ... either in macro or inline
#define LINE_STRING	STRINGIZE(__LINE__)
#define PROFILE_NAME __FILE__ "(" LINE_STRING ")" //FUNCTION_STRING //function not expanding 
#define PROFILE()	Profiler profiler(PROFILE_NAME);
#define PROFILER_DONE()	Profiler::done();
