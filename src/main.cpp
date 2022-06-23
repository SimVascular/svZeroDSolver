#include <iostream>
#include <string>

#include "dofhandler.hpp"
#include "block.hpp"

int main()
{
    DOFHandler dofhandler; // Create an object of MyClass

    // Access attributes and set values
    int size = dofhandler.size();
    int id = dofhandler.register_variable("Hallo1");
    int did = dofhandler.register_variable("Hallo2");
    int eid = dofhandler.register_equation();
    int pid = dofhandler.register_equation();

    // Print attribute values
    std::cout << size << std::endl;
    std::cout << id << std::endl;
    std::cout << dofhandler.size() << std::endl;
    std::cout << "l = { ";
    for (std::string n : dofhandler.variables)
    {
        std::cout << n << ", ";
    }
    std::cout << "};\n";
    std::cout << pid << std::endl;
    std::cout << did << std::endl;

    Block::Parameters params;
    Block block = Block(params, "BlockNameHallo");
    block.setup_dofs(dofhandler);

    return 0;
}