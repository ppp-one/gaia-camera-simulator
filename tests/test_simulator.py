import pytest
import time
import logging
from gaia_camera_simulator import GaiaCameraSimulator

# Set up logging to see the output during tests
logging.basicConfig(level=logging.INFO)


class TestGaiaCameraSimulator:
    @pytest.fixture(autouse=True)
    def setup_and_teardown(self):
        """Fixture to create a GaiaCameraSimulator instance and ensure cleanup."""
        self.simulator = GaiaCameraSimulator(observatory_config_name="default")
        if hasattr(self, "simulator") and self.simulator is not None:
            self.simulator.cleanup()

    def test_start_simulator(self):
        """Test starting the simulator."""
        self.simulator.start()
        assert self.simulator.is_running() is True
        time.sleep(1)  # Give some time for the process to start
        assert self.simulator.is_running() is True

    def test_stop_simulator(self):
        """Test stopping the simulator."""
        self.simulator.start()
        assert self.simulator.is_running() is True
        self.simulator.stop()
        assert self.simulator.is_running() is False

    def test_double_start(self):
        """Test starting the simulator twice."""
        self.simulator.start()
        assert self.simulator.is_running() is True
        self.simulator.start()
        assert self.simulator.is_running() is True
        self.simulator.stop()

    def test_double_stop(self):
        """Test stopping the simulator twice."""
        self.simulator.start()
        assert self.simulator.is_running() is True
        self.simulator.stop()
        assert self.simulator.is_running() is False
        self.simulator.stop()
        assert self.simulator.is_running() is False

    def test_context_manager(self):
        """Test the context manager functionality."""
        with GaiaCameraSimulator(observatory_config_name="default") as sim:
            sim.start()
            assert sim.is_running() is True
        assert sim.is_running() is False  # Should be stopped after exiting context

    def test_cleanup_on_exit(self):
        """Test cleanup on exit."""
        self.simulator.start()
        assert self.simulator.is_running() is True
        del self.simulator  # Force deletion
        time.sleep(1)  # Give some time for cleanup
