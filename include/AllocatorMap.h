#pragma once

#include "Common/Exception.h"
#include <string>
#include <map>
#include <memory>

template<typename Type>
struct Allocator {
	virtual std::shared_ptr<Type> operator()() {
		throw Common::Exception() << "not implemented";
	}
};

template<typename Type, typename ChildType>
struct AllocatorType : public Allocator<Type> {
	virtual std::shared_ptr<Type> operator()() {
		return std::dynamic_pointer_cast<Type>(std::make_shared<ChildType>());
	}
};

template<typename Type>
struct AllocatorMap {
	typedef std::map<std::string, std::shared_ptr<Allocator<Type>>> Map;
	Map map;
	
	std::shared_ptr<Type> operator()(const std::string& name);

	template<typename ChildType>
	void add(std::string name) {
		map[name] = std::make_shared<AllocatorType<Type, ChildType>>();
	}
};

template<typename Type>
std::shared_ptr<Type> AllocatorMap<Type>::operator()(const std::string& name) {
	for (typename Map::value_type &value : map) {
		if (value.first == name) return (*value.second)();
	}
	throw Common::Exception() << "unknown " << name;
}

