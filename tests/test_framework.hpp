/**
 * SeisProc Test Framework
 * 
 * A lightweight unit testing framework for SeisProc components.
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <cmath>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <map>

namespace realdetect {
namespace test {

// Test result
struct TestResult {
    std::string name;
    bool passed;
    std::string message;
    double duration_ms;
};

// Test case
struct TestCase {
    std::string name;
    std::function<void()> func;
};

// Global test registry
class TestRegistry {
public:
    static TestRegistry& instance() {
        static TestRegistry reg;
        return reg;
    }
    
    void addTest(const std::string& suite, const std::string& name, 
                 std::function<void()> func) {
        suites_[suite].push_back({name, func});
    }
    
    std::vector<TestResult> runAll() {
        std::vector<TestResult> results;
        
        for (const auto& [suite_name, tests] : suites_) {
            std::cout << "\n=== " << suite_name << " ===" << std::endl;
            
            for (const auto& test : tests) {
                TestResult result;
                result.name = suite_name + "::" + test.name;
                
                auto start = std::chrono::high_resolution_clock::now();
                
                try {
                    current_test_ = &result;
                    test.func();
                    result.passed = !failed_;
                    if (failed_) {
                        result.message = failure_message_;
                    }
                } catch (const std::exception& e) {
                    result.passed = false;
                    result.message = std::string("Exception: ") + e.what();
                } catch (...) {
                    result.passed = false;
                    result.message = "Unknown exception";
                }
                
                auto end = std::chrono::high_resolution_clock::now();
                result.duration_ms = std::chrono::duration<double, std::milli>(
                    end - start).count();
                
                // Reset state
                failed_ = false;
                failure_message_.clear();
                
                // Print result
                if (result.passed) {
                    std::cout << "  ✓ " << test.name 
                              << " (" << std::fixed << std::setprecision(2) 
                              << result.duration_ms << "ms)" << std::endl;
                } else {
                    std::cout << "  ✗ " << test.name << ": " 
                              << result.message << std::endl;
                }
                
                results.push_back(result);
            }
        }
        
        return results;
    }
    
    std::vector<TestResult> runSuite(const std::string& suite_name) {
        std::vector<TestResult> results;
        
        auto it = suites_.find(suite_name);
        if (it == suites_.end()) {
            std::cerr << "Suite not found: " << suite_name << std::endl;
            return results;
        }
        
        std::cout << "\n=== " << suite_name << " ===" << std::endl;
        
        for (const auto& test : it->second) {
            TestResult result;
            result.name = suite_name + "::" + test.name;
            
            auto start = std::chrono::high_resolution_clock::now();
            
            try {
                current_test_ = &result;
                test.func();
                result.passed = !failed_;
                if (failed_) {
                    result.message = failure_message_;
                }
            } catch (const std::exception& e) {
                result.passed = false;
                result.message = std::string("Exception: ") + e.what();
            }
            
            auto end = std::chrono::high_resolution_clock::now();
            result.duration_ms = std::chrono::duration<double, std::milli>(
                end - start).count();
            
            failed_ = false;
            failure_message_.clear();
            
            if (result.passed) {
                std::cout << "  ✓ " << test.name << std::endl;
            } else {
                std::cout << "  ✗ " << test.name << ": " << result.message << std::endl;
            }
            
            results.push_back(result);
        }
        
        return results;
    }
    
    void fail(const std::string& message) {
        failed_ = true;
        failure_message_ = message;
    }
    
    bool isFailed() const { return failed_; }

private:
    TestRegistry() = default;
    std::map<std::string, std::vector<TestCase>> suites_;
    TestResult* current_test_ = nullptr;
    bool failed_ = false;
    std::string failure_message_;
};

// Test registration helper
struct TestRegistrar {
    TestRegistrar(const std::string& suite, const std::string& name,
                  std::function<void()> func) {
        TestRegistry::instance().addTest(suite, name, func);
    }
};

// Assertion macros
#define ASSERT_TRUE(cond) \
    if (!(cond)) { \
        std::ostringstream oss; \
        oss << "ASSERT_TRUE failed: " #cond " at " << __FILE__ << ":" << __LINE__; \
        realdetect::test::TestRegistry::instance().fail(oss.str()); \
        return; \
    }

#define ASSERT_FALSE(cond) \
    if (cond) { \
        std::ostringstream oss; \
        oss << "ASSERT_FALSE failed: " #cond " at " << __FILE__ << ":" << __LINE__; \
        realdetect::test::TestRegistry::instance().fail(oss.str()); \
        return; \
    }

#define ASSERT_EQ(a, b) \
    if ((a) != (b)) { \
        std::ostringstream oss; \
        oss << "ASSERT_EQ failed: " << (a) << " != " << (b) \
            << " at " << __FILE__ << ":" << __LINE__; \
        realdetect::test::TestRegistry::instance().fail(oss.str()); \
        return; \
    }

#define ASSERT_NE(a, b) \
    if ((a) == (b)) { \
        std::ostringstream oss; \
        oss << "ASSERT_NE failed: " << (a) << " == " << (b) \
            << " at " << __FILE__ << ":" << __LINE__; \
        realdetect::test::TestRegistry::instance().fail(oss.str()); \
        return; \
    }

#define ASSERT_LT(a, b) \
    if (!((a) < (b))) { \
        std::ostringstream oss; \
        oss << "ASSERT_LT failed: " << (a) << " >= " << (b) \
            << " at " << __FILE__ << ":" << __LINE__; \
        realdetect::test::TestRegistry::instance().fail(oss.str()); \
        return; \
    }

#define ASSERT_LE(a, b) \
    if (!((a) <= (b))) { \
        std::ostringstream oss; \
        oss << "ASSERT_LE failed: " << (a) << " > " << (b) \
            << " at " << __FILE__ << ":" << __LINE__; \
        realdetect::test::TestRegistry::instance().fail(oss.str()); \
        return; \
    }

#define ASSERT_GT(a, b) \
    if (!((a) > (b))) { \
        std::ostringstream oss; \
        oss << "ASSERT_GT failed: " << (a) << " <= " << (b) \
            << " at " << __FILE__ << ":" << __LINE__; \
        realdetect::test::TestRegistry::instance().fail(oss.str()); \
        return; \
    }

#define ASSERT_GE(a, b) \
    if (!((a) >= (b))) { \
        std::ostringstream oss; \
        oss << "ASSERT_GE failed: " << (a) << " < " << (b) \
            << " at " << __FILE__ << ":" << __LINE__; \
        realdetect::test::TestRegistry::instance().fail(oss.str()); \
        return; \
    }

#define ASSERT_NEAR(a, b, eps) \
    if (std::abs((a) - (b)) > (eps)) { \
        std::ostringstream oss; \
        oss << "ASSERT_NEAR failed: |" << (a) << " - " << (b) << "| = " \
            << std::abs((a) - (b)) << " > " << (eps) \
            << " at " << __FILE__ << ":" << __LINE__; \
        realdetect::test::TestRegistry::instance().fail(oss.str()); \
        return; \
    }

#define ASSERT_THROW(expr, exc_type) \
    { bool caught = false; \
      try { expr; } catch (const exc_type&) { caught = true; } catch (...) {} \
      if (!caught) { \
          std::ostringstream oss; \
          oss << "ASSERT_THROW failed: expected " #exc_type \
              << " at " << __FILE__ << ":" << __LINE__; \
          realdetect::test::TestRegistry::instance().fail(oss.str()); \
          return; \
      } \
    }

#define ASSERT_NO_THROW(expr) \
    try { expr; } catch (const std::exception& e) { \
        std::ostringstream oss; \
        oss << "ASSERT_NO_THROW failed: " << e.what() \
            << " at " << __FILE__ << ":" << __LINE__; \
        realdetect::test::TestRegistry::instance().fail(oss.str()); \
        return; \
    }

// Test definition macro
#define TEST(suite, name) \
    void test_##suite##_##name(); \
    static realdetect::test::TestRegistrar registrar_##suite##_##name( \
        #suite, #name, test_##suite##_##name); \
    void test_##suite##_##name()

// Print summary
inline void printSummary(const std::vector<TestResult>& results) {
    int passed = 0, failed = 0;
    double total_time = 0;
    
    for (const auto& r : results) {
        if (r.passed) passed++;
        else failed++;
        total_time += r.duration_ms;
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test Summary" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Passed: " << passed << "/" << (passed + failed) << std::endl;
    std::cout << "Failed: " << failed << std::endl;
    std::cout << "Total time: " << std::fixed << std::setprecision(2) 
              << total_time << "ms" << std::endl;
    
    if (failed > 0) {
        std::cout << "\nFailed tests:" << std::endl;
        for (const auto& r : results) {
            if (!r.passed) {
                std::cout << "  - " << r.name << ": " << r.message << std::endl;
            }
        }
    }
    
    std::cout << "========================================" << std::endl;
}

} // namespace test
} // namespace realdetect
