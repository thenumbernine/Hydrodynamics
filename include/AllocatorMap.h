#pragma once

#include "Common/Exception.h"
#include <string>
#include <map>
#include <functional>

template<typename Type>
struct AllocatorMap {
	typedef std::map<std::string, std::function<Type *()>> Map;
	Map map;
	
	Type *create(const std::string &name);
};

template<typename Type>
Type *AllocatorMap<Type>::create(const std::string &name) {
	for (const typename Map::value_type &value : map) {
		if (value.first == name) return value.second();
	}
	throw Common::Exception() << "unknown " << name;
}

