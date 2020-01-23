function simulationMethod = ResolveIncrementMethod(simulationMethodString, numberOfCoefficients)

switch simulationMethodString
    case "gouincrementconttankov"
        simulationMethod = GOUIncrementContTankov();
    case "gouincrementcufarosabino"
        simulationMethod = GOUIncrementCufaroSabino();
    case "gouincrementcufarosabinorejection"
        simulationMethod = GOUIncrementCufaroSabinoRejection(numberOfCoefficients);
    case "gouincrementqdz"
        simulationMethod = GOUIncrementQDZ();
    case "bgousymmetricincrementconttankov"
        simulationMethod = BGOUSymmetricIncrementContTankov();
    case "bgousymmetricincrementcufarosabino"
        simulationMethod = BGOUSymmetricIncrementCufaroSabino();
    case "bgousymmetricincrementcufarosabinorejection"
        simulationMethod = BGOUSymmetricIncrementCufaroSabinoRejection(numberOfCoefficients);
    case "bgousymmetricincrementqdz"
        simulationMethod = BGOUSymmetricIncrementQDZ();
    otherwise
        error("Increment method not implemented");
end
