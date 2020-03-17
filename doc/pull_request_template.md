<!--
This template should be included in all pull requests. Items in the list should
either be completed by the original author or explicitly dismissed by one of the
IBAMR principal developers.

IBAMR is a community effort and it wouldn't exist without people contributing
code. Thanks in advance for helping to make IBAMR better!
-->

### IBAMR Pull Request Checklist
- [ ] Does the test suite pass?
- [ ] Was clang-format run?
- [ ] Were relevant issues cited? Please remember to add `Fixes #12345` to close
      the issue automatically if we are fixing the problem.
- [ ] Is this a change others will want to know about? If so, then has a
      changelog entry been added?
- [ ] If new data structures have been added to a class then have they been set
      up appropriately for restarts? If so, ensure that the restart version
      number is incremented.
- [ ] Does this change include a bug fix or new functionality? If so, a new test
      or tests should be added. New tests should run quickly (less than a minute
      in release mode). If possible, an older test should gain a new option so
      that we do not need to compile more test executables.
- [ ] Did you (if your account has permission to do so) set relevant labels on
      GitHub for the pull request?
