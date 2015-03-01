#ifndef ARGUMENTHANDLER_H
#define ARGUMENTHANDLER_H

#include <string>
#include <vector>
#include <map>

enum ArgumentFlag {REQUIRED_ARGUMENT, OPTIONAL_ARGUMENT, SWITCH, FLEXIBLE_SWITCH};

template <typename T>
class ArgumentWrapper {
public:
	ArgumentWrapper(T& handle) : handle(&handle) { }
	void SetValue(T value) { *handle = value; }
	ArgumentFlag flag;
	T defaultValue;
	T switchValue;
	std::string description;
private:
	T* handle;
};

class ArgumentInterface {
public:
	virtual ~ArgumentInterface() {}
	virtual std::string GetDescription(unsigned int index) = 0;
	virtual ArgumentFlag GetFlag(unsigned int index) = 0;
	virtual void SetToDefaultValue(unsigned int index) = 0;
	virtual void SetToSwitchValue(unsigned int index) = 0;
	virtual void SetToPassedString(unsigned int index, std::string value) = 0;
};

template <typename T>
class ArgumentList : public ArgumentInterface {
public:
	unsigned int Size() const { return arglist.size(); }
	void Add(ArgumentWrapper<T> arg);
	
	std::string GetDescription(unsigned int index);
	ArgumentFlag GetFlag(unsigned int index);
	void SetToDefaultValue(unsigned int index);
	void SetToSwitchValue(unsigned int index);
	void SetToPassedString(unsigned int index, std::string value);
	void SetToPassedValue(unsigned int index, T value);
private:
	std::vector <ArgumentWrapper<T> > arglist;
};

class ArgumentMap {
private:
	struct MapIndex {
		MapIndex() { }
		MapIndex(char mapchar, int index) : mapchar(mapchar), index(index) { }
		char mapchar;
		int index;
	};
	std::map <std::string, MapIndex> argMap;
	std::vector <std::string> argOrder;
	
	void AddMapping(std::string name, char mapchar, int index);
	ArgumentInterface& Selector(char typechar);
	std::string ArgumentDescriptionBuilder(ArgumentInterface& argList, unsigned int index, std::string name, char typechar);

	ArgumentList <double> argListDouble;
	ArgumentList <int> argListInt;
	ArgumentList <std::string> argListString;
	ArgumentList <bool> argListBool;
	ArgumentList <unsigned> argListUnsigned;
public:
	void PassArgumentValue(std::string arg, std::string value);
	void PassSwitch(std::string arg);
	void PassNothing(std::string arg);
	std::string GetDescription(std::string arg);
	
	void Add(std::string name, ArgumentWrapper<double> variable);
	void Add(std::string name, ArgumentWrapper<int> variable);
	void Add(std::string name, ArgumentWrapper<std::string> variable);
	void Add(std::string name, ArgumentWrapper<bool> variable);
	void Add(std::string name, ArgumentWrapper<unsigned> variable);
};

class ArgumentHandler {
private:
	struct BoolString {
		BoolString() : set(false) {}
		BoolString(std::string str) : strval(str), set(true) {}
		BoolString& operator= (std::string str) { strval = str; set = true; return *this; }
		BoolString& operator= (bool flag) { strval = ""; set = flag; return *this; }
		operator std::string() const { return strval; }
		operator bool() const { return set; }
	private:
		std::string strval;
		bool set;
	};
	ArgumentMap expectedArgs;
	std::vector <std::string> argNames;
	std::string programName;
	std::map <std::string, BoolString> readArgs;

	void ParsingReader(int argc, char* argv[]);
public:
	void Parse(int argc, char* argv[]);
	void PrintUsageSyntax();

	template <typename T> void AddRequiredArgument(std::string name, T& variable, std::string description);
	template <typename T, typename U> void AddOptionalArgument(std::string name, T& variable, U defaultValue, std::string description);
	template <typename T, typename U> void AddFlexibleSwitch(std::string name, T& variable, U switchOnValue, U defaultValue, std::string description);
	void AddSwitch(std::string name, bool& variable, std::string description);
};

template <typename T>
void ArgumentHandler::AddRequiredArgument(std::string name, T& variable, std::string description) {
	ArgumentWrapper<T> temp(variable);
	temp.flag = REQUIRED_ARGUMENT;
	temp.description = description;
	expectedArgs.Add(name, temp);
	argNames.push_back(name);
}

template <typename T, typename U>
void ArgumentHandler::AddOptionalArgument(std::string name, T& variable, U defaultValue, std::string description) {
	ArgumentWrapper<T> temp(variable);
	temp.flag = OPTIONAL_ARGUMENT;
	temp.defaultValue = defaultValue;
	temp.description = description;
	expectedArgs.Add(name, temp);
	argNames.push_back(name);
}

template <typename T, typename U>
void ArgumentHandler::AddFlexibleSwitch(std::string name, T& variable, U switchOnValue, U defaultValue, std::string description) {
	ArgumentWrapper<T> temp(variable);
	temp.flag = FLEXIBLE_SWITCH;
	temp.defaultValue = defaultValue;
	temp.switchValue = switchOnValue;
	temp.description = description;
	expectedArgs.Add(name, temp);
	argNames.push_back(name);
}

#endif
