#pragma once

#include "Common/Exception.h"

#include <math.h>
#include <sys/time.h>

#include <iostream>
#include <string>
#include <map>

struct Profiler {
	static double getTime() {
		timeval t;
		gettimeofday(&t, NULL);
		return (double)t.tv_sec + (double)t.tv_usec * 1e-6;
	}

	struct ProfileEntry {
		double duration;
		double durationExcludingChildren;
		const char *name;
		const char *parentName;
		ProfileEntry()
		: duration(HUGE_VAL)
		, durationExcludingChildren(HUGE_VAL)
		, name(NULL)
		, parentName(NULL)
		{}

		std::map<const char *, bool> childNames;

		void print(int spaces = 0) {
			for (int i = 0; i < spaces; ++i) std::cout << " ";
			std::cout 
				<< duration 
				<< " (" << durationExcludingChildren << ")   " 
				<< name << std::endl;
			++spaces;	
			std::for_each(childNames.begin(), childNames.end(), [&](std::map<const char *, bool>::value_type v) {
				profileMap[v.first].print(spaces);
			});
		}
	};

	static std::map<const char *, ProfileEntry> profileMap;
	
	const char *name;
	double startTime, endTime, duration, subDurations;

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
		endTime = getTime();
		duration = endTime - startTime;
		if (lastProfiler) {
			lastProfiler->subDurations += duration;
		}

		rootName = name;
		currentProfiler = lastProfiler;

		ProfileEntry &entry = profileMap[name];
		entry.name = name;
		entry.duration = std::min(entry.duration, duration);
		entry.durationExcludingChildren = std::min(entry.durationExcludingChildren, duration - subDurations);
		if (lastProfiler) {
			entry.parentName = lastProfiler->name;
			profileMap[lastProfiler->name].childNames[entry.name] = true;
		}
	}

	static void done() {
		if (!rootName) throw Exception() << "need to run Profiler at least once";
		profileMap[rootName].print();
	}
};


//http://stackoverflow.com/questions/3859944/combining-string-literals-and-integer-constants
#define STRINGIZE_DETAIL_(v)	#v
#define STRINGIZE(v) STRINGIZE_DETAIL_(v)
#define FUNCTION_STRING	STRINGIZE(__func__)	//not working ... either in macro or inline
#define LINE_STRING	STRINGIZE(__LINE__)
#define PROFILE_NAME __FILE__ "(" LINE_STRING "): " FUNCTION_STRING 
#define PROFILE()	Profiler profiler(PROFILE_NAME);

