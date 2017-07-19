#include <Model/tests/test_resource_model.h>
#include <gtest/gtest.h>
#include "Simulation/simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/schedule_section.h"

using namespace ::Simulation::SimulatorInterfaces::DriverFileWriters::DriverParts::ECLDriverParts;

namespace {

class DriverPartScheduleTest : public ::testing::Test, public TestResources::TestResourceModel {
protected:
    DriverPartScheduleTest(){
        QStringList *driver_file_contents = Utilities::FileHandling::ReadFileToStringList(settings_simulator_->driver_file_path());
        schedule_ = new Schedule(model_->wells(), settings_model_->control_times(),driver_file_contents);
    }
    virtual ~DriverPartScheduleTest(){}

    Schedule *schedule_;
};

TEST_F(DriverPartScheduleTest, Constructor) {
    //std::cout << schedule_->GetPartString().toStdString() << std::endl;
}

}
