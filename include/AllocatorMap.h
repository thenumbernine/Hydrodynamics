#pragma once

#include <string>
#include <map>
#include <memory>
#include <functional>

template<typename Type>
using AllocatorMap = std::map<std::string, std::function<std::shared_ptr<Type>()>>;
